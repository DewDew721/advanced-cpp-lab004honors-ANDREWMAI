#include "project/DynamicProgramming.hpp"

#include <stdexcept>
#include <climits>
#include <algorithm>

std::string buildParenthesization(
    const std::vector<std::vector<int>>& split,
    int i,
    int j) {
    if (i == j) {
        return "A" + std::to_string(i);
    }

    int k = split[i][j];
    std::string left = buildParenthesization(split, i, k);
    std::string right = buildParenthesization(split, k + 1, j);
    return "(" + left + " x " + right + ")";
}

MatrixChainResult matrixChainOrder(
    const std::vector<int>& dimensions) {
    if (dimensions.size() < 2) {
        throw std::invalid_argument("At least two dimensions are required for one matrix chain");
    }

    for (int dim : dimensions) {
        if (dim <= 0) {
            throw std::invalid_argument("All matrix dimensions must be positive");
        }
    }

    int n = static_cast<int>(dimensions.size()) - 1;
    std::vector<std::vector<long long>> cost(n + 1, std::vector<long long>(n + 1, 0));
    std::vector<std::vector<int>> split(n + 1, std::vector<int>(n + 1, -1));

    for (int length = 2; length <= n; ++length) {
        for (int i = 1; i <= n - length + 1; ++i) {
            int j = i + length - 1;
            cost[i][j] = LLONG_MAX;

            for (int k = i; k < j; ++k) {
                long long q = cost[i][k] + cost[k + 1][j]
                    + static_cast<long long>(dimensions[i - 1]) * dimensions[k] * dimensions[j];

                if (q < cost[i][j]) {
                    cost[i][j] = q;
                    split[i][j] = k;
                }
            }
        }
    }

    return {cost[1][n], split, buildParenthesization(split, 1, n)};
}

AlignmentResult globalSequenceAlignment(
    const std::string& first,
    const std::string& second,
    int matchScore,
    int mismatchScore,
    int gapScore) {
    size_t m = first.size();
    size_t n = second.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    for (size_t i = 0; i <= m; ++i) {
        dp[i][0] = static_cast<int>(i) * gapScore;
    }
    for (size_t j = 0; j <= n; ++j) {
        dp[0][j] = static_cast<int>(j) * gapScore;
    }

    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int matchOrMismatch = (first[i - 1] == second[j - 1]) ? matchScore : mismatchScore;
            dp[i][j] = std::max({
                dp[i - 1][j - 1] + matchOrMismatch,
                dp[i - 1][j] + gapScore,
                dp[i][j - 1] + gapScore
            });
        }
    }

    std::string alignedFirst, alignedSecond;
    size_t i = m;
    size_t j = n;

    while (i > 0 && j > 0) {
        int matchOrMismatch = (first[i - 1] == second[j - 1]) ? matchScore : mismatchScore;
        int diagonal = dp[i - 1][j - 1] + matchOrMismatch;
        int up = dp[i - 1][j] + gapScore;
        int left = dp[i][j - 1] + gapScore;

        if (up >= diagonal && up >= left) {
            alignedFirst = first[i - 1] + alignedFirst;
            alignedSecond = '-' + alignedSecond;
            --i;
        } else if (left >= diagonal && left >= up) {
            alignedFirst = '-' + alignedFirst;
            alignedSecond = second[j - 1] + alignedSecond;
            --j;
        } else {
            alignedFirst = first[i - 1] + alignedFirst;
            alignedSecond = second[j - 1] + alignedSecond;
            --i;
            --j;
        }
    }

    while (i > 0) {
        alignedFirst = first[i - 1] + alignedFirst;
        alignedSecond = '-' + alignedSecond;
        --i;
    }
    while (j > 0) {
        alignedFirst = '-' + alignedFirst;
        alignedSecond = second[j - 1] + alignedSecond;
        --j;
    }

    return {dp[m][n], alignedFirst, alignedSecond};
}
