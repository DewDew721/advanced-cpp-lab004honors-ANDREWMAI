#include "project/DynamicProgramming.hpp"

#include <stdexcept>

std::string buildParenthesization(
    const std::vector<std::vector<int>>& split,
    int i,
    int j) {

    // TODO: Implement recursive reconstruction of the optimal
    // matrix-chain parenthesization using the split table.
    // Base case: if i == j, return "A" + std::to_string(i).
    // Recursive case: split at split[i][j] and combine both subchains.

    if (i == j)
    {
        return "A" + std::to_string(i);
    }
    else
    {
        int k = split[i][j];
        std::string left = buildParenthesization(split, i, k);
        std::string right = buildParenthesization(split, k + 1, j);
        return "(" + left + " x " + right + ")";
    }
}

MatrixChainResult matrixChainOrder(
    const std::vector<int>& dimensions) {

    // TODO: Implement Matrix Chain-Product using bottom-up dynamic programming.
    // dimensions = {p0, p1, ..., pn} represents n matrices.
    // Create a cost table and a split table, evaluate increasing chain lengths,
    // store the best split, and reconstruct the optimal parenthesization.
    // Validate that at least one matrix exists and all dimensions are positive.
    // Do not hard-code the sample answer.

    int n = dimensions.size() - 1;
    std::vector<std::vector<long long>> cost(n + 1, std::vector<long long>(n + 1, 0));
    std::vector<std::vector<int>> split(n + 1, std::vector<int>(n + 1, 0));

    for (int length = 2; length <= n; ++length) {
        for (int i = 1; i <= n - length + 1; ++i) {
            int j = i + length - 1;
            cost[i][j] = LLONG_MAX;
            for (int k = i; k < j; ++k) {
                long long q = cost[i][k] + cost[k + 1][j] + static_cast<long long>(dimensions[i - 1]) * dimensions[k] * dimensions[j];
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

    // TODO: Implement global sequence alignment using dynamic programming.
    // Build and initialize the DP table, evaluate diagonal/up/left transitions,
    // then backtrack from the final cell to reconstruct both aligned strings.
    // Use '-' for gaps and support arbitrary scoring parameters.
    // Do not hard-code the DNA/text examples.

    size_t m = first.size();
    size_t n = second.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    // Initialize the DP table
    for (size_t i = 0; i <= m; ++i) {
        dp[i][0] = i * gapScore;
    }
    for (size_t j = 0; j <= n; ++j) {
        dp[0][j] = j * gapScore;
    }

    // Fill the DP table
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int match = dp[i - 1][j - 1] + (first[i - 1] == second[j - 1] ? matchScore : mismatchScore);
            int deleteOp = dp[i - 1][j] + gapScore;
            int insertOp = dp[i][j - 1] + gapScore;
            dp[i][j] = std::max({match, deleteOp, insertOp});
        }
    }

    // Backtrack to reconstruct the aligned strings
    std::string alignedFirst, alignedSecond;
    size_t i = m, j = n;
    while (i > 0 && j > 0) {
        if (first[i - 1] == second[j - 1]) {
            alignedFirst = first[i - 1] + alignedFirst;
            alignedSecond = second[j - 1] + alignedSecond;
            --i;
            --j;
        } else if (dp[i - 1][j] > dp[i][j - 1]) {
            alignedFirst = first[i - 1] + alignedFirst;
            alignedSecond = '-' + alignedSecond;
            --i;
        } else {
            alignedFirst = '-' + alignedFirst;
            alignedSecond = second[j - 1] + alignedSecond;
            --j;
        }
    }

    // Handle remaining characters
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
