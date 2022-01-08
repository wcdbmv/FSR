#include <algorithm/max_diff.hpp>

#include <numeric>


namespace {

template <typename Container, typename Callback>
Float transform_max(const Container& a, const Container& b, Callback callback) {
	return std::transform_reduce(
		a.begin(),
		a.end(),
		b.begin(),
		Float{},
		[](Float acc, Float cur) { return std::max(acc, cur); },
		callback
	);
}

} // namespace


Float max_diff(const Vector& a, const Vector& b) {
	return transform_max(a, b, [](Float ai, Float bi) { return std::abs((ai - bi) / ai); });
}

Float max_diff(const Matrix& a, const Matrix& b) {
	return transform_max(a, b, [](const auto& ai, const auto& bi) { return max_diff(ai, bi); });
}
