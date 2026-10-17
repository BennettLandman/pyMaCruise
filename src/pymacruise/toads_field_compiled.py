"""Strict source-order compiled kernels for TOADS polynomial field fitting."""
from __future__ import annotations
import numpy as np
from numba import njit


@njit(cache=True, fastmath=False)
def _basis(point, dimensions, degree):
    x = np.float32(point[0]) / np.float32(dimensions[0] - 1)
    y = np.float32(point[1]) / np.float32(dimensions[1] - 1)
    z = np.float32(point[2]) / np.float32(dimensions[2] - 1)
    terms = np.zeros(20, dtype=np.float32)
    terms[0] = np.float32(1.0)
    terms[1:4] = (x, y, z)
    if degree >= 2:
        terms[4] = np.float32(np.float32(2.0) * x * x - np.float32(1.0))
        terms[5] = np.float32(np.float32(2.0) * y * y - np.float32(1.0))
        terms[6] = np.float32(np.float32(2.0) * z * z - np.float32(1.0))
        terms[7] = np.float32(terms[1] * terms[2])
        terms[8] = np.float32(terms[2] * terms[3])
        terms[9] = np.float32(terms[3] * terms[1])
    if degree >= 3:
        terms[10] = np.float32(np.float32(np.float32(4.0) * x * x * x) - np.float32(3.0) * x)
        terms[11] = np.float32(np.float32(np.float32(4.0) * y * y * y) - np.float32(3.0) * y)
        terms[12] = np.float32(np.float32(np.float32(4.0) * z * z * z) - np.float32(3.0) * z)
        terms[13] = np.float32(terms[4] * terms[2])
        terms[14] = np.float32(terms[1] * terms[5])
        terms[15] = np.float32(terms[5] * terms[3])
        terms[16] = np.float32(terms[2] * terms[6])
        terms[17] = np.float32(terms[6] * terms[1])
        terms[18] = np.float32(terms[3] * terms[4])
        terms[19] = np.float32(np.float32(terms[1] * terms[2]) * terms[3])
    return terms[: (4, 10, 20)[degree - 1]].copy()


@njit(cache=True, fastmath=False)
def accumulate_field(image, mem, mask, centers, field_values, degree, background, optimized_classes):
    dimensions = image.shape
    n_parameters = (4, 10, 20)[degree - 1]
    rhs = np.zeros(n_parameters, dtype=np.float64)
    normal = np.zeros((n_parameters, n_parameters), dtype=np.float64)
    nx, ny, nz = dimensions
    for x in range(0, nx, 3):
        for y in range(0, ny, 3):
            for z in range(0, nz, 3):
                if not mask[x, y, z]:
                    continue
                basis = _basis((x, y, z), dimensions, degree)
                for klass in range(mem.shape[-1]):
                    if klass == background:
                        continue
                    corrected = np.float32(field_values[x, y, z] * image[x, y, z])
                    residual = np.float32(corrected - centers[klass])
                    distance = np.float32(residual * residual)
                    optimized_weight = np.float32(1.0)
                    if distance > np.float32(1e-30) and optimized_classes[klass]:
                        # PowerTable exponent=1, range [0,2], step=1e-6.
                        # Preserve its double quantization before float32 division.
                        powered = float(distance)
                        if powered <= 2.0:
                            step = np.float64(np.float32(1e-6))
                            powered = np.trunc(powered / step + 0.5) * step
                        optimized_weight = np.float32(powered / float(distance))

                    membership_squared = np.float32(mem[x, y, z, klass] * mem[x, y, z, klass])
                    vector_term = np.float32(membership_squared * optimized_weight)
                    vector_term = np.float32(vector_term * image[x, y, z])
                    vector_term = np.float32(vector_term * centers[klass])
                    matrix_term = np.float32(membership_squared * optimized_weight)
                    matrix_term = np.float32(matrix_term * image[x, y, z])
                    matrix_term = np.float32(matrix_term * image[x, y, z])
                    for row in range(n_parameters):
                        rhs[row] += float(np.float32(vector_term * basis[row]))
                        row_term = np.float32(matrix_term * basis[row])
                        for column in range(n_parameters):
                            normal[row, column] += float(np.float32(row_term * basis[column]))

    return normal, rhs


@njit(cache=True, fastmath=False)
def evaluate_field(field_values, mask, coefficients, degree):
    dimensions = field_values.shape
    nx, ny, nz = dimensions
    n_parameters = (4, 10, 20)[degree - 1]
    mean = np.float32(0.0)
    count = np.float32(0.0)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if not mask[x, y, z]:
                    continue
                basis = _basis((x, y, z), dimensions, degree)
                value = np.float32(0.0)
                for index in range(n_parameters):
                    value = np.float32(float(coefficients[index]) * float(basis[index]) + float(value))
                field_values[x, y, z] = value
                mean = np.float32(mean + value)
                count = np.float32(count + np.float32(1.0))
    mean = np.float32(mean / count)
    offset = np.float32(np.float32(1.0) - mean)
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if mask[x, y, z]:
                    field_values[x, y, z] = np.float32(field_values[x, y, z] + offset)

