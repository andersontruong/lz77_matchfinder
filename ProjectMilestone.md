# Introduction: what is the problem, why is it interesting?
This project focuses on improving match finding on SCL's LZ77 implementation. Currently, the existing implementations either have no window (`LZ77Encoder`) or a sliding window (`LZ77SlidingWindowEncoder`) to bound memory usage. The newer `LZ77SlidingWindowEncoder` relies on a `MatchFinderBase` to produce an LZ77 sequence at a given position.

The current match finder (`HashBasedMatchFinder`) uses a hash table where each entry contains a hash chain. The hash chain size is limited to a fixed value (`max_chain_length`). This bounds the number of similar matches in order to reduce computation time. The hash keys are currently a hash of the first minimum number of bytes of a sequence (e.g. 3). If we encounter more sequences than the `max_chain_length`, we remove older entries. This leads to less optimal matches when we encounter a large number of repeats. For example, given `hash_length = 2` and `max_chain_length = 64`, if we encounter a stream of "AB...AB..." with more than 64 repeats of sequences starting with "AB", we will lose the first few sequences that may provide a better (longer) match. We could just increase `max_chain_length`, but this will inevitably lead to more computation.

For `HashBasedMatchFinder`, the hash table is typically sparse. On repeated or similar sequences, we end up having a significant number of hash collisions, which leads to the previous problem. To improve upon this match finder, we need a method to reduce the number of hash collisions to better balance the table and reduce sparsity.

# Literature/Code review: summarize any existing papers/implementations that you referred to
- Morphing Match Chain https://fastcompression.blogspot.com/p/mmc-morphing-match-chain.html
- Progressive Hash Series (same principles as MMC) https://fastcompression.blogspot.com/2011/10/progressive-hash-series-new-method-to.html

Yann Collet has blog posts on two different hash-chain based match finders that improve upon the current SCL implementation. [Morphing Match Chain (MMC)](https://fastcompression.blogspot.com/p/mmc-morphing-match-chain.html) is one that hierarchically stores positions based on increasing hash lengths without repeats. This technically lazily builds a suffix tree, where at each level N characters have already been matches.

As we perform searches on a given position, we keep checking for longer matches and promote previous matches. For example, let's look at the input sequence "AABAABCABABCDABCABCD". If we are searching for "ABC" at position 13, first we can look the hash table for a length N=1 match of "A". Then we want to find positions for matches of (N+1) = 2, or "AB. We find positions 7, 9 in the "A" chain that match to "AB". Now we add position 13 to the "A" chain and remove 7, 9 from the "A" chain. Now we go to the "AB" chain, add the previously removed positions 7, 9 to "AB" if they aren't already there (7 was already there), then finally add 13 to the "AB" chain as the newest one. Now we look at the "AB" chain for matches of (N+1) = 3, or "ABC". We find positions 4, 9 that match to "ABC" (4 was already on the "AB" chain from before). We add 13 to the "AB" chain and remove 4, 9. Now we go to the "ABC" chain, add 4, 9, and finally add 13. Now we look at the "ABC" chain for (N+1) = 4, or "ABCA". We see that position 4 matches to "ABCA" and we remove it from "ABC". Then we go to "ABCA" chain, add 4, then 13. We look at "ABCA" chain for "ABCAB", find and remove 4 again. We go to "ABCAB", add 4 and add 13. We look at "ABCAB" chain for "ABCABC" and don't find anything. We add 13 to the new "ABCABC" chain. Then, using the last valid match (position 4 for "ABCAB") we use that match (pos 4, length 5). Then using the rest of the characters of the match at 13 "BCAB", we add them to their respective single character chains 14->"B", 15->"C", 16->"A", 17->"B".

On the previous example, given min_hash_length = 1, the standard `HashBasedMatchFinder` has a maximum hash chain length of 7 (for chains "A" and "B"), average hash chain length of 5, and uses 4 hash chains. Using MMC, the longest hash chain is 5 (for chain "B" since we don't really start searches at "B"), average hash chain length is 1.94, and we use 17 hash chains. For future searches, MMC will be much faster since we can simply descend the hierarchy and not waste time on redundant early matches. With much shorter average chain lengths, we are less likely to run into the maximum chain length that may remove more optimal older matches. This effectively reduces the hash collisions and better utilizes the table.

# Methods: what do you plan to implement as part of this project? What end result do you expect to achieve and how will you evaluate it qualitatively and quantitatively?
In this project, I plan on subclassing `HashBasedMatchFinder` and implementing improved `find_best_match_at_position` and `find_best_match` functions that allow for much larger hash chain lengths to improve compression without greatly increasing processing time. I expect much faster compression speeds (even in Python) with the same window size constraints, particularly on repeating data.

We may also need better metrics, such as number of lookups and hash table utilization. I will work on exposing these counts from the `MatchFinderBase` class.

# Progress report: what have you already finished (please include code link where relevant)? What is the plan for the remaining weeks?
Currently I have worked through MMC on a sample input to verify what the final hash table should look like. I am using this to aid my implementation. I have also set up the basic structure of the new `MorphingMatchChainBasedMatchFinder` by inheriting `HashBasedMatchFinder`. My plan for the rest of the quarter is to finish implementing this and verify the compression efficiency against a range of `max_chain_length` on the existing benchmarks in `lz77_sliding_window.py`.

Lazy MMC
  - Current LZ77 Sliding Window is asymptotically optimal as window size goes to infinity [Wyner-Ziv, 94]
  - For current testing, we can just massively increase the window size and measure latency

Pseudo-Suffix Trie
- May be cycles due to hash collisions (apparent in BootstrapJS benchmark)

TODO:
1. Fix the tunneling problem of promotion (multi-level promotion)
2. Instead of just adding the single matches at the end, do secondary promotions
3. Add lazy lookahead
4. Improve runtime?
