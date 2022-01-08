#pragma once

#include <vector>

template <typename T>
using Container = std::vector<T>;

using Float = float;
using Size = Container<Float>::size_type;

using Vector = Container<Float>;

using Matrix = Container<Vector>;
