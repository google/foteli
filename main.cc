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

#include <OpenImageIO/imageio.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>

#include "foteli.h"

int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input.png> <output.png>\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const char* output_filename = argv[2];

  const auto start_time_load = std::chrono::high_resolution_clock::now();
  auto input = OIIO::ImageInput::open(input_filename);
  if (!input) {
    fprintf(stderr, "Failed to open \"%s\".\n", input_filename);
    return EXIT_FAILURE;
  }
  const OIIO::ImageSpec& spec = input->spec();
  if (spec.nchannels != 3) {
    fprintf(stderr, "Expected an image with 3 channels, got %d\n",
            spec.nchannels);
    return EXIT_FAILURE;
  }
  std::vector<std::unique_ptr<float[]>> channels;
  for (int c = 0; c < spec.nchannels; ++c) {
    std::unique_ptr<float[]> channel_data(new float[spec.width * spec.height]);
    if (!input->read_image(/*subimage=*/0, /*miplevel=*/0, c, c + 1,
                           OIIO::TypeDesc::FLOAT, channel_data.get())) {
      fprintf(stderr, "Error reading \"%s\": %s\n", input_filename,
              input->geterror().c_str());
      return EXIT_FAILURE;
    }
    channels.push_back(std::move(channel_data));
  }
  input->close();
  const auto end_time_load = std::chrono::high_resolution_clock::now();

  fprintf(
      stderr, "Loaded a %dx%d, %d-channel image in %ldms.\n", spec.width,
      spec.height, spec.nchannels,
      static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(
                            end_time_load - start_time_load)
                            .count()));

  const auto start_time_compute = std::chrono::high_resolution_clock::now();
  FoteliParams* const params = FoteliCreateParams();
  FoteliSetFloatBuffers(params, channels[0].get(), channels[1].get(),
                        channels[2].get());
  FoteliSetImageSize(params, spec.width, spec.height);
  if (const char* error_message = FoteliValidateParams(params)) {
    fprintf(stderr, "Foteli error: %s\n", error_message);
    return EXIT_FAILURE;
  }
  FoteliToneMap(params);
  FoteliDestroyParams(params);
  const auto end_time_compute = std::chrono::high_resolution_clock::now();

  fprintf(
      stderr, "Computed in %ldms.\n",
      static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(
                            end_time_compute - start_time_compute)
                            .count()));

  const auto start_time_save = std::chrono::high_resolution_clock::now();
  auto output = OIIO::ImageOutput::create(output_filename);
  if (!output) {
    fprintf(stderr, "Failed to create output for \"%s\"\n", output_filename);
    return EXIT_FAILURE;
  }
  OIIO::ImageSpec output_spec(spec.width, spec.height, spec.nchannels,
                              OIIO::TypeDesc::UINT16);
  if (!output->open(output_filename, output_spec)) {
    fprintf(stderr,
            "Failed to open \"%s\" to write a 16-bit %dx%dx%d image: %s\n",
            output_filename, spec.width, spec.height, spec.nchannels,
            output->geterror().c_str());
    return EXIT_FAILURE;
  }
  std::unique_ptr<uint16_t[]> row(new uint16_t[spec.width * spec.nchannels]);
  for (int y = 0; y < spec.height; ++y) {
    for (int x = 0; x < spec.width; ++x) {
      for (int c = 0; c < spec.nchannels; ++c) {
        row[x * spec.nchannels + c] = static_cast<uint16_t>(
            .5f + 65535.f * channels[c][y * spec.width + x]);
      }
    }
    if (!output->write_scanline(y, 0, OIIO::TypeDesc::UINT16, row.get())) {
      fprintf(stderr, "Failed to write row %d: %s\n", y,
              output->geterror().c_str());
      return EXIT_FAILURE;
    }
  }
  if (!output->close()) {
    fprintf(stderr, "Failed to flush \"%s\": %s\n", output_filename,
            output->geterror().c_str());
    return EXIT_FAILURE;
  }
  const auto end_time_save = std::chrono::high_resolution_clock::now();
  fprintf(
      stderr, "Saved in %ldms.\n",
      static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(
                            end_time_save - start_time_save)
                            .count()));
}
