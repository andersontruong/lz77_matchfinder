import argparse
from scl.compressors.lz77_sliding_window import MatchFinderBase, LZ77SlidingWindowEncoder

class MorphingMatchChainBasedMatchFinder(HashBasedMatchFinder):
    def __init__(
        self,
        hash_length=4,
        hash_table_size=1000000,
        max_chain_length=64,
        lazy=True,
        minimum_match_length=4,
    ):
        """
        Initializes the morphing match chain based match finder
        """
        super().__init__(
            hash_length=hash_length,
            hash_table_size=hash_table_size,
            max_chain_length=max_chain_length,
            lazy=lazy,
            minimum_match_length=minimum_match_length
        )

    def add_to_hashtable(self, position, data):
        """Updates the hash table with the new position of a byte sequence. Maintains chain length limit."""
        if len(data) >= self.hash_length:
            h = self._hash(data[: self.hash_length])
            if len(self.hash_table[h]) == self.max_chain_length:
                self.hash_table[h].pop(0)
            self.hash_table[h].append(position)
            self.next_position_to_hash = position + 1

    def find_best_match(self, lookahead_buffer):
        """
        Find the best match for the current position in the lookahead buffer
        within the sliding window using the hash table.

        Parameters:
        - lookahead_buffer (bytearray): The buffer containing the data that we're trying to find a match for.

        Returns:
        - tuple: A tuple containing the number of literals, the match position, and the match length.
        """
        i = 0

        # we can look up to the end of the lookahead buffer minus the end part
        # that we can't yet compute the hash for
        num_positions_to_consider = len(lookahead_buffer) - self.hash_length + 1
        if num_positions_to_consider <= 0:
            # this means lookahead_buffer is smaller than the hash length
            # so we just return the length of the buffer as literals
            return (len(lookahead_buffer), 0, 0)

        for i in range(num_positions_to_consider):
            # Ensure all bytes are hashed and added to the hash table
            # We want to index until the end of the window and then also
            # the lookahead buffer until i
            for j in range(max(self.next_position_to_hash, self.window.start), self.window.end + i):
                print(f'\t\tAdding pos {j} to hash')
                # create the key
                data = bytearray()
                for k in range(j, j + self.hash_length):
                    data.append(self.window.get_byte_window_plus_lookahead(k, lookahead_buffer))
                self.add_to_hashtable(j, data)

            print(f'\n\t{i}')

        raise NotImplementedError


if __name__ == "__main__":
    # Provide a simple CLI interface below for convenient experimentation
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file", required=True, type=str)
    parser.add_argument("-o", "--output", help="output file", required=True, type=str)
    parser.add_argument(
        "-w", "--window_init", help="initialize window from file (like zstd dictionary)", type=str
    )

    # constants
    BLOCKSIZE = 100_000  # encode in 100 KB blocks

    args = parser.parse_args()

    initial_window = None
    if args.window_init is not None:
        with open(args.window_init, "rb") as f:
            initial_window = list(f.read())

    match_finder = MorphingMatchChainBasedMatchFinder(hash_length=1, minimum_match_length=1)
    encoder = LZ77SlidingWindowEncoder(match_finder, initial_window=initial_window)
    encoder.encode_file(args.input, args.output, block_size=BLOCKSIZE)
