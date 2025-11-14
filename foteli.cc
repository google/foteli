// Copyright 2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "foteli.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "foteli.cc"
#include <hwy/foreach_target.h>
#include <hwy/highway.h>

#include <hwy/contrib/algo/transform-inl.h>
#include <hwy/contrib/math/math-inl.h>
#include <hwy/contrib/thread_pool/thread_pool.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#ifndef FOTELI_DEFINE_FOTELI_PARAMS_
#define FOTELI_DEFINE_FOTELI_PARAMS_

struct FoteliParams {
  float* red = nullptr;
  float* green = nullptr;
  float* blue = nullptr;
  int width = 0, height = 0;
  int row_stride = 0;
};

#endif

namespace {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

constexpr int kDownsampling = 128;
constexpr float kDefaultIntensityTarget = 255;
// Base-2 log.
constexpr float kLogDefaultIntensityTarget = 7.994353436858858;
constexpr float kLog10000 = 13.287712379549449;

struct ImageParams {
  explicit ImageParams(const FoteliParams& params)
      : width(params.width),
        height(params.height),
        row_stride(params.row_stride) {}
  ImageParams(const int width, const int height)
      : width(width), height(height), row_stride(width) {}
  int width, height, row_stride;
};

enum Access { kReadOnly, kReadWrite };

template <Access>
struct Image {};

template <>
struct Image<kReadOnly> {
  const float* pixels;
  ImageParams params;
  const float* operator[](const int y) const {
    return pixels + y * params.row_stride;
  }
};

template <>
struct Image<kReadWrite> {
  float* pixels;
  ImageParams params;
  const float* operator[](const int y) const {
    return pixels + y * params.row_stride;
  }
  float* operator[](const int y) { return pixels + y * params.row_stride; }

  operator Image<kReadOnly>() const {
    return {.pixels = pixels, .params = params};
  }
};

struct OwnedImage {
  static OwnedImage CopyFrom(const Image<kReadOnly>& other) {
    OwnedImage result(other.params.width, other.params.height);
    for (int y = 0; y < other.params.height; ++y) {
      std::copy_n(other[y], other.params.width, result.image[y]);
    }
    return result;
  }
  explicit OwnedImage(const int width, const int height)
      : pixels(new float[width * height]),
        image{.pixels = pixels.get(), .params = ImageParams(width, height)} {}
  std::unique_ptr<float[]> pixels;
  Image<kReadWrite> image;
  float* operator[](const int y) { return image[y]; }
  const float* operator[](const int y) const { return image[y]; }
};

template <typename D, typename V>
HWY_ATTR V PQToLinear(D d, const V pq) {
  const V x = hn::MulAdd(pq, pq, pq);
  static constexpr float kP[5] = {5500.34862f, 26455.3172f, 7386.02301f,
                                  -62.3553089f, 2.62975656f};
  static constexpr float kQ[5] = {2.67718770f, -33.9078883f, 174.364667f,
                                  -428.736818f, 421.350107f};
  V yp = hn::Set(d, kP[0]);
  V yq = hn::Set(d, kQ[0]);
  for (int i = 1; i < 5; ++i) {
    yp = hn::MulAdd(yp, x, hn::Set(d, kP[i]));
    yq = hn::MulAdd(yq, x, hn::Set(d, kQ[i]));
  }
  return hn::Div(yp, yq);
}

template <typename D, typename V>
HWY_ATTR V LinearToSRGB(D d, const V linear) {
  static constexpr float kP[5] = {7.352629620e-01f, 1.474205315e+00f,
                                  3.903842876e-01f, 5.287254571e-03f,
                                  -5.135152395e-04f};
  static constexpr float kQ[5] = {2.424867759e-02f, 9.258482155e-01f,
                                  1.340816930e+00f, 3.036675394e-01f,
                                  1.004519624e-02f};

  const V if_small = hn::Mul(hn::Set(d, 12.92f), linear);

  const V x = hn::Sqrt(linear);
  V yp = hn::Set(d, kP[0]);
  V yq = hn::Set(d, kQ[0]);
  for (int i = 1; i < 5; ++i) {
    yp = hn::MulAdd(yp, x, hn::Set(d, kP[i]));
    yq = hn::MulAdd(yq, x, hn::Set(d, kQ[i]));
  }
  const V if_large = hn::Div(yp, yq);

  return hn::IfThenElse(hn::Le(linear, hn::Set(d, 0.04045f / 12.92f)), if_small,
                        if_large);
}

template <typename D, typename V>
HWY_ATTR V ComputeLuminance(D d, const V r, const V g, const V b) {
  return hn::Max(hn::Set(d, 1e-12f),
                 hn::MulAdd(hn::Set(d, 0.2627f), r,
                            hn::MulAdd(hn::Set(d, 0.6780f), g,
                                       hn::Mul(hn::Set(d, 0.0593f), b))));
}

int DivCeil(const int a, const int b) { return (a + b - 1) / b; }

HWY_ATTR OwnedImage DownsampledLuminances(const Image<kReadOnly>& red,
                                          const Image<kReadOnly>& green,
                                          const Image<kReadOnly>& blue,
                                          const int width, const int height,
                                          hwy::ThreadPool& pool) {
  const int downsampled_width = DivCeil(width, kDownsampling);
  const int downsampled_height = DivCeil(height, kDownsampling);
  OwnedImage result(downsampled_width, downsampled_height);

  HWY_FULL(float) d;
  pool.Run(
      0, downsampled_height, [&](const int y, const int /* worker */) HWY_ATTR {
        for (int x = 0; x < downsampled_width; ++x) {
          auto max_luminance = hn::Set(d, .5f * kDefaultIntensityTarget);
          for (int ky = 0; ky < kDownsampling; ++ky) {
            if (y * kDownsampling + ky >= height) break;
            const float* r = red[y * kDownsampling + ky] + x * kDownsampling;
            const float* g = green[y * kDownsampling + ky] + x * kDownsampling;
            const float* b = blue[y * kDownsampling + ky] + x * kDownsampling;
            const int max_kx =
                std::min(kDownsampling, width - x * kDownsampling);
            int kx;
            for (kx = 0; kx + Lanes(d) <= max_kx; kx += Lanes(d)) {
              auto red = hn::LoadU(d, &r[kx]);
              auto green = hn::LoadU(d, &g[kx]);
              auto blue = hn::LoadU(d, &b[kx]);
              max_luminance =
                  hn::Max(max_luminance, ComputeLuminance(d, red, green, blue));
            }
            if (kx <= max_kx) {
              const int n = max_kx - kx;
              const auto no = hn::Zero(d);
              auto red = hn::LoadNOr(no, d, &r[kx], n);
              auto green = hn::LoadNOr(no, d, &g[kx], n);
              auto blue = hn::LoadNOr(no, d, &b[kx], n);
              max_luminance =
                  hn::Max(max_luminance, ComputeLuminance(d, red, green, blue));
            }
          }
          result[y][x] = hn::ReduceMax(d, max_luminance);
        }
        hn::Transform(d, result[y], downsampled_width,
                      [](auto d, auto max_luminance)
                          HWY_ATTR { return hn::Log2(d, max_luminance); });
      });

  return result;
}

HWY_ATTR void Blur(Image<kReadWrite> image) {
  HWY_FULL(float) d;
  static constexpr std::array<float, 3> kWeights = {.375f, .25f, .0625f};
  OwnedImage blurred_horizontally(image.params.width, image.params.height);
  for (int y = 0; y < image.params.height; ++y) {
    for (int x = 0; x < image.params.width; ++x) {
      blurred_horizontally[y][x] = 0;
      for (int dx = -2; dx <= 2; ++dx) {
        const int clamped_x =
            std::max(0, std::min(image.params.width - 1, x + dx));
        blurred_horizontally[y][x] +=
            image[y][clamped_x] * kWeights[std::abs(dx)];
      }
    }
  }
  for (int y = 0; y < image.params.height; ++y) {
    hn::Generate(d, image[y], image.params.width,
                 [](auto d, auto index) HWY_ATTR { return hn::Zero(d); });
    for (int dy = -2; dy <= 2; ++dy) {
      const auto weight = hn::Set(d, kWeights[std::abs(dy)]);
      const int clamped_y =
          std::max(0, std::min(image.params.height - 1, y + dy));
      hn::Transform1(
          d, image[y], image.params.width, blurred_horizontally[clamped_y],
          [&weight](auto d, auto image, auto blurred_horizontally) HWY_ATTR {
            return hn::MulAdd(blurred_horizontally, weight, image);
          });
    }
  }
}

HWY_ATTR OwnedImage Upsample(const Image<kReadOnly>& image,
                             hwy::ThreadPool& pool) {
  HWY_FULL(float) d;
  using V = decltype(hn::Zero(d));
  OwnedImage upsampled_horizontally(2 * image.params.width,
                                    image.params.height);
  const auto BoundX = [&image](int x) {
    return std::max(0, std::min(image.params.width - 1, x));
  };
  pool.Run(0, image.params.height, [&](const int y, const int /* worker */) {
    const float* const in_row = image[y];
    float* const out_row = upsampled_horizontally[y];

    for (int x = 0; x < image.params.width; ++x) {
      out_row[2 * x] = in_row[x];
      out_row[2 * x + 1] =
          0.5625f * (in_row[x] + in_row[BoundX(x + 1)]) -
          0.0625f * (in_row[BoundX(x - 1)] + in_row[BoundX(x + 2)]);
    }
  });

  OwnedImage upsampled(2 * image.params.width, 2 * image.params.height);
  const auto BoundY = [&image](int y) {
    return std::max(0, std::min(image.params.height - 1, y));
  };
  pool.Run(
      0, image.params.height,
      [&](const int y, const int /* worker */) HWY_ATTR {
        const float* const in_rows[4] = {
            upsampled_horizontally[BoundY(y - 1)],
            upsampled_horizontally[y],
            upsampled_horizontally[BoundY(y + 1)],
            upsampled_horizontally[BoundY(y + 2)],
        };
        float* HWY_RESTRICT const out_rows[2] = {
            upsampled[2 * y],
            upsampled[2 * y + 1],
        };

        int x;
        for (x = 0; x + Lanes(d) <= upsampled_horizontally.image.params.width;
             x += Lanes(d)) {
          const V in[] = {
              hn::LoadU(d, &in_rows[0][x]),
              hn::LoadU(d, &in_rows[1][x]),
              hn::LoadU(d, &in_rows[2][x]),
              hn::LoadU(d, &in_rows[3][x]),
          };
          hn::StoreU(in[1], d, &out_rows[0][x]);
          hn::StoreU(
              hn::MulAdd(hn::Set(d, 0.5625f), hn::Add(in[1], in[2]),
                         hn::Mul(hn::Set(d, -0.0625f), hn::Add(in[0], in[3]))),
              d, &out_rows[1][x]);
        }
        if (x < upsampled_horizontally.image.params.width) {
          const int n = upsampled_horizontally.image.params.width - x;
          V in[] = {
              hn::LoadN(d, &in_rows[0][x], n),
              hn::LoadN(d, &in_rows[1][x], n),
              hn::LoadN(d, &in_rows[2][x], n),
              hn::LoadN(d, &in_rows[3][x], n),
          };
          hn::StoreN(in[1], d, &out_rows[0][x], n);
          hn::StoreN(
              hn::MulAdd(hn::Set(d, 0.5625f), hn::Add(in[1], in[2]),
                         hn::Mul(hn::Set(d, -0.0625f), hn::Add(in[0], in[3]))),
              d, &out_rows[1][x], n);
        }
      });
  return upsampled;
}

HWY_ATTR void ApplyToneMapping(Image<kReadWrite> r, Image<kReadWrite> g,
                               Image<kReadWrite> b,
                               const Image<kReadOnly>& blurred_luminances,
                               hwy::ThreadPool& pool) {
  HWY_FULL(float) d;

  auto ComputeRatio = [&d](auto blurred_luminances, auto r, auto g,
                           auto b) HWY_ATTR {
    const auto log_local_max = blurred_luminances + hn::Set(d, 1);
    const auto luminance = ComputeLuminance(d, r, g, b);
    const auto log_luminance = hn::Min(log_local_max, hn::Log2(d, luminance));
    const auto log_knee = hn::Mul(
        hn::Set(d, kLogDefaultIntensityTarget),
        hn::MulAdd(hn::Set(d, -0.85f),
                   hn::Div(hn::Sub(log_local_max,
                                   hn::Set(d, kLogDefaultIntensityTarget)),
                           hn::Set(d, kLog10000 - kLogDefaultIntensityTarget)),
                   hn::Set(d, 1)));
    const auto second_segment_position = hn::Div(
        hn::Sub(log_luminance, log_knee), hn::Sub(log_local_max, log_knee));
    const auto log_new_luminance = IfThenElse(
        hn::Lt(log_luminance, log_knee), log_luminance,
        hn::MulAdd(
            second_segment_position,
            hn::MulAdd(
                hn::Sub(hn::Set(d, kLogDefaultIntensityTarget), log_knee),
                second_segment_position, hn::Sub(log_knee, log_luminance)),
            log_luminance));
    const auto new_luminance = hn::Exp2(d, log_new_luminance);
    return hn::Div(new_luminance,
                   hn::Mul(luminance, hn::Set(d, kDefaultIntensityTarget)));
  };

  pool.Run(
      0, r.params.height, [&](const int y, const int /* worker */) HWY_ATTR {
        float* HWY_RESTRICT row_r = r[y];
        float* HWY_RESTRICT row_g = g[y];
        float* HWY_RESTRICT row_b = b[y];
        int x;
        for (x = 0; x + Lanes(d) <= r.params.width; x += Lanes(d)) {
          const auto blurred_lums = hn::LoadU(d, &blurred_luminances[y][x]);
          const auto r = hn::LoadU(d, &row_r[x]);
          const auto g = hn::LoadU(d, &row_g[x]);
          const auto b = hn::LoadU(d, &row_b[x]);
          const auto ratio = ComputeRatio(blurred_lums, r, g, b);
          hn::StoreU(r * ratio, d, &row_r[x]);
          hn::StoreU(g * ratio, d, &row_g[x]);
          hn::StoreU(b * ratio, d, &row_b[x]);
        }
        if (x < r.params.width) {
          const int n = r.params.width - x;
          const auto blurred_lums = hn::LoadN(d, &blurred_luminances[y][x], n);
          const auto r = hn::LoadN(d, &row_r[x], n);
          const auto g = hn::LoadN(d, &row_g[x], n);
          const auto b = hn::LoadN(d, &row_b[x], n);
          const auto ratio = ComputeRatio(blurred_lums, r, g, b);
          hn::StoreN(r * ratio, d, &row_r[x], n);
          hn::StoreN(g * ratio, d, &row_g[x], n);
          hn::StoreN(b * ratio, d, &row_b[x], n);
        }
      });
}

template <typename D, typename V>
HWY_ATTR void Rec2020To709(D d, V& r, V& g, V& b) {
  const V new_r = hn::MulAdd(
      hn::Set(d, 1.6605f), r,
      hn::MulAdd(hn::Set(d, -0.5876f), g, hn::Mul(hn::Set(d, -0.0728f), b)));
  const V new_g = hn::MulAdd(
      hn::Set(d, -0.1246f), r,
      hn::MulAdd(hn::Set(d, 1.1329f), g, hn::Mul(hn::Set(d, -0.0084f), b)));
  b = hn::MulAdd(
      hn::Set(d, -0.0181f), r,
      hn::MulAdd(hn::Set(d, -0.1006f), g, hn::Mul(hn::Set(d, 1.1187f), b)));
  r = new_r;
  g = new_g;
}

template <typename D, typename V>
HWY_ATTR void GamutMap(D d, V& r, V& g, V& b,
                       const std::array<float, 3>& primaries_luminances,
                       const float preserve_saturation = .4f) {
  const V luminance =
      hn::MulAdd(hn::Set(d, primaries_luminances[0]), r,
                 hn::MulAdd(hn::Set(d, primaries_luminances[1]), g,
                            hn::Mul(hn::Set(d, primaries_luminances[2]), b)));

  // Desaturate out-of-gamut pixels. This is done by mixing each pixel
  // with just enough gray of the target luminance to make all
  // components non-negative.
  // - For saturation preservation, if a component is still larger than
  // 1 then the pixel is normalized to have a maximum component of 1.
  // That will reduce its luminance.
  // - For luminance preservation, getting all components below 1 is
  // done by mixing in yet more gray. That will desaturate it further.
  V gray_mix_saturation = hn::Zero(d);
  V gray_mix_luminance = hn::Zero(d);
  for (const V& val : {r, g, b}) {
    const V val_minus_gray = hn::Sub(val, luminance);
    const V inv_val_minus_gray = hn::Div(
        hn::Set(d, 1.f), hn::IfThenElse(hn::Eq(val_minus_gray, hn::Zero(d)),
                                        hn::Set(d, 1), val_minus_gray));
    const V val_over_val_minus_gray = hn::Mul(val, inv_val_minus_gray);
    gray_mix_saturation =
        hn::IfThenElse(hn::Ge(val_minus_gray, hn::Zero(d)), gray_mix_saturation,
                       hn::Max(gray_mix_saturation, val_over_val_minus_gray));
    gray_mix_luminance = hn::Max(
        gray_mix_luminance,
        hn::IfThenElse(hn::Le(val_minus_gray, hn::Zero(d)), gray_mix_saturation,
                       hn::Sub(val_over_val_minus_gray, inv_val_minus_gray)));
  }
  const V gray_mix =
      hn::Clamp(hn::MulAdd(hn::Set(d, preserve_saturation),
                           hn::Sub(gray_mix_saturation, gray_mix_luminance),
                           gray_mix_luminance),
                hn::Zero(d), hn::Set(d, 1.f));
  for (V* const ch : {&r, &g, &b}) {
    *ch = hn::MulAdd(gray_mix, hn::Sub(luminance, *ch), *ch);
  }
  const V max_clr = hn::Max(hn::Max(hn::Set(d, 1.f), r), hn::Max(g, b));
  const V normalizer = hn::Div(hn::Set(d, 1.f), max_clr);
  for (V* const ch : {&r, &g, &b}) {
    *ch = hn::Mul(*ch, normalizer);
  }
}

HWY_ATTR void FoteliToneMapImpl(const FoteliParams* const params) {
  HWY_FULL(float) d;
  hwy::ThreadPool pool(std::min<size_t>(hwy::ThreadPool::MaxThreads(), 16));

  Image<kReadWrite> red = {.pixels = params->red,
                           .params = ImageParams(*params)};
  Image<kReadWrite> green = {.pixels = params->green,
                             .params = ImageParams(*params)};
  Image<kReadWrite> blue = {.pixels = params->blue,
                            .params = ImageParams(*params)};

  for (auto* plane : {&red, &green, &blue}) {
    pool.Run(0, plane->params.height,
             [&](const int y, const int /* worker */) HWY_ATTR {
               hn::Transform(d, (*plane)[y], plane->params.width,
                             [](auto d, auto pixels)
                                 HWY_ATTR { return PQToLinear(d, pixels); });
             });
  }

  auto downsampled_luminances = DownsampledLuminances(
      red, green, blue, params->width, params->height, pool);

  Blur(downsampled_luminances.image);

  OwnedImage blurred_luminances = std::move(downsampled_luminances);
  for (int downsampling = kDownsampling; downsampling > 1; downsampling >>= 1) {
    blurred_luminances = Upsample(blurred_luminances.image, pool);
  }

  ApplyToneMapping(red, green, blue, blurred_luminances.image, pool);

  static constexpr std::array<float, 3> kSRGBLuminances = {.2126f, .7152f,
                                                           .0722f};
  pool.Run(0, params->height,
           [&](const int y, const int /* worker */) HWY_ATTR {
             int x;
             for (x = 0; x + Lanes(d) <= params->width; x += Lanes(d)) {
               auto r = hn::LoadU(d, &red[y][x]);
               auto g = hn::LoadU(d, &green[y][x]);
               auto b = hn::LoadU(d, &blue[y][x]);
               Rec2020To709(d, r, g, b);
               GamutMap(d, r, g, b, kSRGBLuminances);
               hn::StoreU(r, d, &red[y][x]);
               hn::StoreU(g, d, &green[y][x]);
               hn::StoreU(b, d, &blue[y][x]);
             }
             if (x < params->width) {
               const int n = params->width - x;
               auto r = hn::LoadN(d, &red[y][x], n);
               auto g = hn::LoadN(d, &green[y][x], n);
               auto b = hn::LoadN(d, &blue[y][x], n);
               Rec2020To709(d, r, g, b);
               GamutMap(d, r, g, b, kSRGBLuminances);
               hn::StoreN(r, d, &red[y][x], n);
               hn::StoreN(g, d, &green[y][x], n);
               hn::StoreN(b, d, &blue[y][x], n);
             }
           });

  for (auto* plane : {&red, &green, &blue}) {
    pool.Run(0, plane->params.height,
             [&](const int y, const int /* worker */) HWY_ATTR {
               float* row = (*plane)[y];
               hn::Transform(d, row, plane->params.width,
                             [](auto d, auto pixels)
                                 HWY_ATTR { return LinearToSRGB(d, pixels); });
             });
  }
}

}  // namespace HWY_NAMESPACE

#if HWY_ONCE
HWY_EXPORT(FoteliToneMapImpl);
#endif

}  // namespace

#if HWY_ONCE

FoteliParams* FoteliCreateParams() { return new FoteliParams; }

void FoteliSetFloatBuffers(FoteliParams* params, float* red, float* green,
                           float* blue) {
  params->red = red;
  params->green = green;
  params->blue = blue;
}

void FoteliSetImageSize(FoteliParams* params, int width, int height) {
  params->width = width;
  params->height = height;
  if (params->row_stride == 0) {
    params->row_stride = width;
  }
}

void FoteliSetRowStride(FoteliParams* params, int stride) {
  params->row_stride = stride;
}

const char* FoteliValidateParams(const FoteliParams* params) {
  if (params->red == nullptr) {
    return "red buffer is null";
  }
  if (params->green == nullptr) {
    return "green buffer is null";
  }
  if (params->blue == nullptr) {
    return "blue buffer is null";
  }
  if (params->width == 0 || params->height == 0) {
    return "size not specified";
  }
  if (params->row_stride < params->width) {
    return "stride set to less than image width";
  }
  return nullptr;
}

void FoteliDestroyParams(FoteliParams* params) { delete params; }

void FoteliToneMap(const FoteliParams* const params) {
  HWY_DYNAMIC_DISPATCH(FoteliToneMapImpl)(params);
}

#endif
