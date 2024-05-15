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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

struct FoteliParams {
  float* red = nullptr;
  float* green = nullptr;
  float* blue = nullptr;
  int width = 0, height = 0;
  int row_stride = 0;
};

namespace {

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

#pragma omp declare simd
float PQToLinear(const float pq) {
  const float x = pq * pq + pq;
  static constexpr float kP[5] = {5500.34862f, 26455.3172f, 7386.02301f,
                                  -62.3553089f, 2.62975656f};
  static constexpr float kQ[5] = {2.67718770f, -33.9078883f, 174.364667f,
                                  -428.736818f, 421.350107f};
  float yp = kP[0];
  float yq = kQ[0];
  for (int i = 1; i < 5; ++i) {
    yp = yp * x + kP[i];
    yq = yq * x + kQ[i];
  }
  return yp / yq;
}

#pragma omp declare simd
float LinearToSRGB(const float linear) {
  static constexpr float kP[5] = {7.352629620e-01f, 1.474205315e+00f,
                                  3.903842876e-01f, 5.287254571e-03f,
                                  -5.135152395e-04f};
  static constexpr float kQ[5] = {2.424867759e-02f, 9.258482155e-01f,
                                  1.340816930e+00f, 3.036675394e-01f,
                                  1.004519624e-02f};
  if (linear <= 0.04045f / 12.92f) {
    return 12.92f * linear;
  } else {
    const float x = std::sqrt(linear);
    float yp = kP[0];
    float yq = kQ[0];
    for (int i = 1; i < 5; ++i) {
      yp = yp * x + kP[i];
      yq = yq * x + kQ[i];
    }
    return yp / yq;
  }
}

#pragma omp declare simd
float ComputeLuminance(const float r, const float g, const float b) {
  return std::max(1e-12f, 0.2627f * r + 0.6780f * g + 0.0593f * b);
}

int DivCeil(const int a, const int b) { return (a + b - 1) / b; }

OwnedImage DownsampledLuminances(const Image<kReadOnly>& red,
                                 const Image<kReadOnly>& green,
                                 const Image<kReadOnly>& blue, const int width,
                                 const int height) {
  const int downsampled_width = DivCeil(width, kDownsampling);
  const int downsampled_height = DivCeil(height, kDownsampling);
  OwnedImage result(downsampled_width, downsampled_height);

#pragma omp parallel for
  for (int y = 0; y < downsampled_height; ++y) {
    for (int x = 0; x < downsampled_width; ++x) {
      float max_luminance = .5f * kDefaultIntensityTarget;
      for (int ky = 0; ky < kDownsampling; ++ky) {
        if (y * kDownsampling + ky >= height) break;
        const float* r = red[y * kDownsampling + ky] + x * kDownsampling;
        const float* g = green[y * kDownsampling + ky] + x * kDownsampling;
        const float* b = blue[y * kDownsampling + ky] + x * kDownsampling;
        const int max_kx = std::min(kDownsampling, width - x * kDownsampling);
#pragma omp simd reduction(max : max_luminance)
        for (int kx = 0; kx < max_kx; ++kx) {
          max_luminance =
              std::max(max_luminance, ComputeLuminance(r[kx], g[kx], b[kx]));
        }
      }
      result[y][x] = std::log2(max_luminance);
    }
  }

  return result;
}

void Blur(Image<kReadWrite> image) {
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
    for (int x = 0; x < image.params.width; ++x) {
      image[y][x] = 0;
      for (int dy = -2; dy <= 2; ++dy) {
        const int clamped_y =
            std::max(0, std::min(image.params.height - 1, y + dy));
        image[y][x] +=
            blurred_horizontally[clamped_y][x] * kWeights[std::abs(dy)];
      }
    }
  }
}

OwnedImage Upsample(const Image<kReadOnly>& image) {
  OwnedImage upsampled_horizontally(2 * image.params.width,
                                    image.params.height);
  const auto BoundX = [&image](int x) {
    return std::max(0, std::min(image.params.width - 1, x));
  };
#pragma omp parallel for
  for (int y = 0; y < image.params.height; ++y) {
    const float* const in_row = image[y];
    float* const out_row = upsampled_horizontally[y];

    for (int x = 0; x < image.params.width; ++x) {
      out_row[2 * x] = in_row[x];
      out_row[2 * x + 1] =
          0.5625f * (in_row[x] + in_row[BoundX(x + 1)]) -
          0.0625f * (in_row[BoundX(x - 1)] + in_row[BoundX(x + 2)]);
    }
  }

  OwnedImage upsampled(2 * image.params.width, 2 * image.params.height);
  const auto BoundY = [&image](int y) {
    return std::max(0, std::min(image.params.height - 1, y));
  };
#pragma omp parallel for
  for (int y = 0; y < image.params.height; ++y) {
    const float* const in_rows[4] = {
        upsampled_horizontally[BoundY(y - 1)],
        upsampled_horizontally[y],
        upsampled_horizontally[BoundY(y + 1)],
        upsampled_horizontally[BoundY(y + 2)],
    };
    float* const out_rows[2] = {
        upsampled[2 * y],
        upsampled[2 * y + 1],
    };

#pragma omp simd
    for (int x = 0; x < upsampled_horizontally.image.params.width; ++x) {
      out_rows[0][x] = in_rows[1][x];
      out_rows[1][x] = 0.5625f * (in_rows[1][x] + in_rows[2][x]) -
                       0.0625f * (in_rows[0][x] + in_rows[3][x]);
    }
  }
  return upsampled;
}

void ApplyToneMapping(Image<kReadWrite> r, Image<kReadWrite> g,
                      Image<kReadWrite> b,
                      const Image<kReadOnly>& blurred_luminances) {
#pragma omp parallel for
  for (int y = 0; y < r.params.height; ++y) {
    float* row_r = r[y];
    float* row_g = g[y];
    float* row_b = b[y];
#pragma omp simd
    for (int x = 0; x < r.params.width; ++x) {
      const float log_local_max = blurred_luminances[y][x] + 1;
      const float luminance = ComputeLuminance(row_r[x], row_g[x], row_b[x]);
      const float log_luminance = std::min(log_local_max, std::log2(luminance));
      const float log_knee =
          kLogDefaultIntensityTarget *
          (1 - 0.85f * (log_local_max - kLogDefaultIntensityTarget) /
                   (kLog10000 - kLogDefaultIntensityTarget));
      const float second_segment_position =
          (log_luminance - log_knee) / (log_local_max - log_knee);
      const float log_new_luminance =
          log_luminance < log_knee
              ? log_luminance
              : second_segment_position *
                        ((kLogDefaultIntensityTarget - log_knee) *
                             second_segment_position +
                         (log_knee - log_luminance)) +
                    log_luminance;
      const float new_luminance = std::exp2(log_new_luminance);
      const float ratio = new_luminance / (luminance * kDefaultIntensityTarget);

      row_r[x] *= ratio;
      row_g[x] *= ratio;
      row_b[x] *= ratio;
    }
  }
}

#pragma omp declare simd linear(r, g, b)
void Rec2020To709(float* const r, float* const g, float* const b) {
  const float new_r = 1.6605f * *r - 0.5876f * *g - 0.0728f * *b;
  const float new_g = -0.1246f * *r + 1.1329f * *g - 0.0084f * *b;
  *b = -0.0181f * *r - 0.1006f * *g + 1.1187f * *b;
  *r = new_r;
  *g = new_g;
}

#pragma omp declare simd linear(r, g, b) \
    uniform(primaries_luminances, preserve_saturation)
void GamutMap(float* r, float* g, float* b,
              const std::array<float, 3>& primaries_luminances,
              const float preserve_saturation = .4f) {
  const float luminance = primaries_luminances[0] * *r +
                          primaries_luminances[1] * *g +
                          primaries_luminances[2] * *b;

  // Desaturate out-of-gamut pixels. This is done by mixing each pixel
  // with just enough gray of the target luminance to make all
  // components non-negative.
  // - For saturation preservation, if a component is still larger than
  // 1 then the pixel is normalized to have a maximum component of 1.
  // That will reduce its luminance.
  // - For luminance preservation, getting all components below 1 is
  // done by mixing in yet more gray. That will desaturate it further.
  float gray_mix_saturation = 0;
  float gray_mix_luminance = 0;
  for (const float* const ch : {r, g, b}) {
    const float val = *ch;
    const float val_minus_gray = val - luminance;
    const float inv_val_minus_gray =
        1.f / (val_minus_gray == 0 ? 1 : val_minus_gray);
    const float val_over_val_minus_gray = val * inv_val_minus_gray;
    gray_mix_saturation =
        val_minus_gray >= 0
            ? gray_mix_saturation
            : std::max(gray_mix_saturation, val_over_val_minus_gray);
    gray_mix_luminance = std::max(
        gray_mix_luminance, val_minus_gray <= 0
                                ? gray_mix_saturation
                                : val_over_val_minus_gray - inv_val_minus_gray);
  }
  const float gray_mix =
      std::max(0.f, std::min(1.f, preserve_saturation * (gray_mix_saturation -
                                                         gray_mix_luminance) +
                                      gray_mix_luminance));
  for (float* const ch : {r, g, b}) {
    *ch += gray_mix * (luminance - *ch);
  }
  const float max_clr = std::max(std::max(1.f, *r), std::max(*g, *b));
  const float normalizer = 1.f / max_clr;
  for (float* const ch : {r, g, b}) {
    *ch *= normalizer;
  }
}

}  // namespace

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
  Image<kReadWrite> red = {.pixels = params->red,
                           .params = ImageParams(*params)};
  Image<kReadWrite> green = {.pixels = params->green,
                             .params = ImageParams(*params)};
  Image<kReadWrite> blue = {.pixels = params->blue,
                            .params = ImageParams(*params)};

  for (auto* plane : {&red, &green, &blue}) {
#pragma omp parallel for
    for (int y = 0; y < plane->params.height; ++y) {
      float* row = (*plane)[y];
#pragma omp simd
      for (int x = 0; x < plane->params.width; ++x) {
        row[x] = PQToLinear(row[x]);
      }
    }
  }

  auto downsampled_luminances =
      DownsampledLuminances(red, green, blue, params->width, params->height);

  Blur(downsampled_luminances.image);

  OwnedImage blurred_luminances = std::move(downsampled_luminances);
  for (int downsampling = kDownsampling; downsampling > 1; downsampling >>= 1) {
    blurred_luminances = Upsample(blurred_luminances.image);
  }

  ApplyToneMapping(red, green, blue, blurred_luminances.image);

  static constexpr std::array<float, 3> kSRGBLuminances = {.2126f, .7152f,
                                                           .0722f};
#pragma omp parallel for
  for (int y = 0; y < params->height; ++y) {
#pragma omp simd
    for (int x = 0; x < params->width; ++x) {
      Rec2020To709(&red[y][x], &green[y][x], &blue[y][x]);
      GamutMap(&red[y][x], &green[y][x], &blue[y][x], kSRGBLuminances);
    }
  }

  for (auto* plane : {&red, &green, &blue}) {
#pragma omp parallel for
    for (int y = 0; y < plane->params.height; ++y) {
      float* row = (*plane)[y];
#pragma omp simd
      for (int x = 0; x < plane->params.width; ++x) {
        row[x] = LinearToSRGB(row[x]);
      }
    }
  }
}
