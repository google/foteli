/*
   Copyright 2024 Google LLC

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       https://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef FOTELI_H_
#define FOTELI_H_

/**
 * @file foteli.h
 *
 * Main header file for Foteli.
 */

#include "foteli_export.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Opaque structure holding the input parameters for Foteli.
 */
typedef struct FoteliParams FoteliParams;

/**
 * Allocates and initialises an instance of FoteliParams.
 * The result should be freed by @ref FoteliDestroyParams once it is not of any
 * more use, although this is not strictly required if the program is about to
 * terminate and the system will reclaim its memory anyway.
 *
 * @return A pointer to a newly allocated and initialised instance of
 *         FoteliParams.
 */
FOTELI_EXPORT FoteliParams* FoteliCreateParams(void);

/**
 * @param params The parameter struct to update.
 * @param red Pointer to the red plane of the pixel data.
 * @param green Pointer to the green plane of the pixel data.
 * @param blue Pointer to the blue plane of the pixel data.
 */
FOTELI_EXPORT void FoteliSetFloatBuffers(FoteliParams* params, float* red,
                                         float* green, float* blue);

/**
 * Sets the size of the image to tone map, in pixels.
 *
 * @param params The parameter struct to update.
 * @param width Width of the image, in pixels.
 * @param height Height of the image, in pixels.
 */
FOTELI_EXPORT void FoteliSetImageSize(FoteliParams* params, int width,
                                      int height);

/**
 * Optional: sets the stride between rows to something else than its default. By
 * default, the stride is set by @ref FoteliSetImageSize to the width of the
 * image.
 *
 * It is expected that for each plane `p`, pixel data for the pixel at (x, y) is
 * at `p[y * row_stride + x]`.
 *
 * @param params The parameter struct to update.
 * @param row_stride The stride, in pixels (not bytes).
 */
FOTELI_EXPORT void FoteliSetRowStride(FoteliParams* params, int row_stride);

/**
 * Verifies that the parameters are usable. This is expected to be most useful
 * during development. If/once you are sure that you never generate invalid
 * parameter structs, it is fine not to call this.
 *
 * @param params The parameter struct to validate.
 * @return @c NULL if the parameters are valid, otherwise an error message
 *         describing the issue. The value returned by this function points to
 *         static memory and must not be freed by the user.
 */
FOTELI_EXPORT const char* FoteliValidateParams(const FoteliParams* params);

/**
 * Deallocates the given parameter struct.
 *
 * @param params Pointer to the parameter struct that was previously allocated
 *               by @ref FoteliCreateParams.
 */
FOTELI_EXPORT void FoteliDestroyParams(FoteliParams* params);

/**
 * Main Foteli function. The input is expected to be in Rec. 2020 / PQ, in the
 * [0, 1] range. It is modified in-place and becomes sRGB in that same range.
 * (Future releases might make parts of this configurable.)
 *
 * @param params The parameter struct holding all the necessary information.
 */
FOTELI_EXPORT void FoteliToneMap(const FoteliParams* params);

#ifdef __cplusplus
}
#endif

#endif  // FOTELI_H_
