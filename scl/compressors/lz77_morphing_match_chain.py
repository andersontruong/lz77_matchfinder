import argparse
from scl.compressors.lz77_sliding_window import HashBasedMatchFinder, LZ77SlidingWindowEncoder, LZ77SlidingWindowDecoder

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

    def extend_match(
        self, match_start_in_lookahead_buffer, match_pos, lookahead_buffer, left_extension=True
    ):
        """Extend candidate match to the right and left as long as the bytes match.
        Returns the match position and length of the match.

        Extend candidate match to the right, promote matches, prune hash chains

        Parameters:
        lookahead_buffer_pos (int): The index of the candidate start match in the lookahead buffer
        match_pos (int): The position of the candidate start match in the sliding window
        lookahead_buffer (bytearray): The lookahead buffer (this is basically the part not yet
                                      added to window where we are looking for a match).
        left_extension (bool): Whether to try to extend the match to left within the lookahead buffer.

        Returns:
        match_start_in_lookahead_buffer (int): The start index of the best match in the lookahead buffer (might be
                                               different from the candidate start match if we extended left)
        match_pos (int): The position of the best match in the sliding window
        length (int): The length of the best match
        """

        # Right extension:
        # Start comparing bytes from the potential match position in the window and the lookahead buffer
        # until a mismatch is found or we've checked all bytes in the lookahead buffer.
        match_length = 0
        while (
            match_start_in_lookahead_buffer + match_length < len(lookahead_buffer)
            and self.window.get_byte_window_plus_lookahead(
                match_pos + match_length, lookahead_buffer
            )
            == lookahead_buffer[match_start_in_lookahead_buffer + match_length]
        ):
            match_length += 1

        # if left_extension:
        #     # Left extension:
        #     # If two sequences match at [i, i+length] in the lookahead buffer and [match_pos, match_pos+length] in the window,
        #     # we then look at the previous byte (i.e., i-1 in the lookahead buffer and match_pos-1 in the window).
        #     # If they match, we continue moving left until a mismatch is found,
        #     # ensuring we don't try to extend beyond the window's start or the lookahead buffer's start.
        #
        #     # Left extension is useful where we missed a hash match because of a collision in the previous
        #     # step. Note that we only extend left within the lookahead buffer, so we don't tread on the
        #     # already encoded matches in past.
        #
        #     # note that this does not change the offset
        #     while match_pos > self.window.start and (
        #         match_start_in_lookahead_buffer > 0
        #         and lookahead_buffer[match_start_in_lookahead_buffer - 1]
        #         == self.window.get_byte_window_plus_lookahead(match_pos - 1, lookahead_buffer)
        #     ):
        #         match_pos -= 1
        #         match_length += 1
        #         match_start_in_lookahead_buffer -= 1

        return (
            match_start_in_lookahead_buffer,
            match_pos,
            match_length,
        )  # match_start_in_lookahead_buffer is the number of literals

    def add_to_hashtable_keep_next_position(self, position, data):
        """Updates the hash table with the new position of a byte sequence. Maintains chain length limit.
           Supports random insert by not updating the next position."""
        if len(data) >= self.hash_length:
            h = self._hash(data)
            if len(self.hash_table[h]) == self.max_chain_length:
                self.hash_table[h].pop(0)
            self.hash_table[h].append(position)
            # We need to manually update since we insert into the hash table multiple times per match
            # self.next_position_to_hash = position + 1

    def find_best_match_at_position(self, lookahead_buffer, i):
        """
        Find the best match for the current position in the lookahead buffer
        within the sliding window using the hash table.

        Parameters:
        - lookahead_buffer (bytearray): The buffer containing the data that we're trying to find a match for.
        - i (int): The position in the lookahead buffer to find a match for.

        Returns:
        - tuple: A tuple containing the number of literals, the match position, and the match length in the best found match.
                 If no match is found, the match length is 0.
        """
        # MMC 1. Add current position as hash to table
        best_match_pos = 0
        best_length = 0
        best_literals_count = i

        # Retrieve potential match positions from the hash table
        hash_data = bytearray(lookahead_buffer[i : i + self.hash_length])
        # print('Hash data: ', hash_data)
        candidate_positions = self.get_positions_from_hash(hash_data)

        # If candidates exist
        while len(candidate_positions) > 0:
            # print(f'\tHash match to {repr(''.join([chr(x) for x in hash_data]))} - [{self.window.end + i}, {self.window.end + i + self.hash_length + best_length})')
            # print(f'\t\tCandidates: {candidate_positions}')

            # 1. Need to validate if any of them are true matches (collisions)
            # - If so, update best length
            # 2. Get bigger hash and validate against true matches
            # - If they match, promote them

            prev_candidates = candidate_positions

            # N+1 Chain
            # next_hash_data = bytearray(lookahead_buffer[i : i + len(hash_data) + 1])
            if i + len(hash_data) < len(lookahead_buffer):
                next_hash_data = hash_data + bytearray([lookahead_buffer[i + len(hash_data)]])
            next_hash = self._hash(next_hash_data)
            next_candidates = self.get_positions_from_hash(
                next_hash_data
            )

            # If prev_candidates can also be in N+1, remove from N and add to N+1
            # print(f'\tChecking candidates for promotion to {len(next_hash_data)} - {repr(''.join([chr(x) for x in next_hash_data]))} : {prev_candidates}')
            prev_candidates_copy = candidate_positions.copy()
            found_candidate = False
            for candidate in prev_candidates_copy:
                offset = self.window.end + i - candidate
                if offset > self.window.size:
                    print('Larger than window, skipping candidate')
                    continue

                # If we already lazily added it, return literal
                # if self.window.end + i == candidate:
                #     print('Already saw...', i)
                #     return (1, 0, 0)

                # Check if true candidate (to avoid hash collisions)
                candidate_match = bytearray()
                for x in range(len(hash_data)):
                    candidate_match.append(self.window.get_byte_window_plus_lookahead(candidate + x, lookahead_buffer))
                if hash_data == candidate_match:
                    # print('\t\tMatch!!!')
                    best_match_pos = candidate
                    found_candidate = True
                else:
                    # print('\t\tNo match :(')
                    continue

                if candidate + len(hash_data) - self.window.end >= len(lookahead_buffer):
                    break
                candidate_nplus1 = candidate_match
                candidate_nplus1.append(self.window.get_byte_window_plus_lookahead(candidate + len(hash_data), lookahead_buffer))
                # print(next_hash_data, 'vs', candidate_nplus1)
                candidate_hash = self._hash(candidate_nplus1)

                # print(f'\t\tCandidate {candidate} - {repr(''.join([chr(x) for x in candidate_nplus1]))} - {candidate_hash} ?= {next_hash}')
                if candidate_hash == next_hash:
                    # print(f'\t\tPromoting candidate at {candidate}')
                    prev_candidates.remove(candidate)
                    if candidate not in next_candidates:
                        self.add_to_hashtable_keep_next_position(candidate, next_hash_data)

            if found_candidate:
                best_length = len(hash_data)
                print('Adding to length, now', best_length)

            # print(f'\t\tUpdated N={len(hash_data)}) candidates: {self.get_positions_from_hash(hash_data)}')
            # print(f'\t\tUpdated N+1={len(next_hash_data)}) candidates: {self.get_positions_from_hash(next_hash_data)}')

            # Add to current N chain
            # print(f'\tAdded to N chain for {repr(''.join([chr(x) for x in hash_data]))} - [{self.window.end + i}, {self.window.end + i + len(hash_data)})')
            self.add_to_hashtable_keep_next_position(self.next_position_to_hash, hash_data)
            # print(f'\t\tNew chain: {self.get_positions_from_hash(hash_data)}')

            if len(next_hash_data) == len(hash_data):
                # print('Breaking...')
                break

            hash_data = next_hash_data
            candidate_positions = next_candidates

        # Effectively create a new N+1 hash chain
        # print(f'\tCreated new N+1 chain for {repr(''.join([chr(x) for x in hash_data]))} - [{self.window.end + i}, {self.window.end + i + len(hash_data)})')
        self.add_to_hashtable_keep_next_position(self.next_position_to_hash, hash_data)
        self.next_position_to_hash += 1
        # print(f'\t\tNew chain: {self.get_positions_from_hash(hash_data)}')

        return best_literals_count, best_match_pos, best_length

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
            print(f'\n\t{i}')
            for j in range(max(self.next_position_to_hash, self.window.start), self.window.end + i):
                print(f'\t\tAdding pos {j} to hash')
                # create the key
                data = bytearray()
                for k in range(j, j + self.hash_length):
                    data.append(self.window.get_byte_window_plus_lookahead(k, lookahead_buffer))
                self.add_to_hashtable(j, data)

            # MMC: this adds current position to potentially multiple hash chains
            # Also lazily prunes hash chains using new position as gateway
            best_literals_count, best_match_pos, best_length = self.find_best_match_at_position(
                lookahead_buffer, i
            )
            # print(f'\tFound best match with literal count {best_literals_count}, starting {best_match_pos} of length {best_length}')
            data = []
            for k in range(best_match_pos, best_match_pos + best_length):
                data.append(self.window.get_byte_window_plus_lookahead(k, lookahead_buffer))
            # print(f'\tStr: {repr("".join([chr(x) for x in data]))}')
            if best_length >= self.minimum_match_length:
                # TODO: Ignore lazy right now
                return (best_literals_count, best_match_pos, best_length)
                if not self.lazy:
                    return (best_literals_count, best_match_pos, best_length)
                else:
                    # we now keep going forward until the match length keeps increasing
                    # or we reach the end of the lookahead buffer
                    for j in range(i + 1, min(i+4, num_positions_to_consider)):
                        print(f'\tLazy at {j}')
                        (
                            new_best_literals_count,
                            new_best_match_pos,
                            new_best_length,
                        ) = self.find_best_match_at_position(lookahead_buffer, j)
                        if new_best_length > best_length:
                            best_literals_count = new_best_literals_count
                            best_match_pos = new_best_match_pos
                            best_length = new_best_length
                            print(f'\t\t Lazy found best match with literal count {best_literals_count}, starting {best_match_pos} of length {best_length}')
                        else:
                            break
                    print(f'\t\tFinal at {best_match_pos} of length {best_length}, count {best_literals_count}')
                    return (best_literals_count, best_match_pos, best_length)

        return (num_positions_to_consider, 0, 0)


if __name__ == "__main__":
    # Provide a simple CLI interface below for convenient experimentation
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--decompress", help="decompress", action="store_true")
    parser.add_argument("-i", "--input", help="input file", required=True, type=str)
    parser.add_argument("-o", "--output", help="output file", required=True, type=str)
    parser.add_argument(
        "-w", "--window_init", help="initialize window from file (like zstd dictionary)", type=str
    )
    parser.add_argument(
        "-c", "--max_chain_length", help="set max chain length", type=int
    )

    # constants
    BLOCKSIZE = 100_000  # encode in 100 KB blocks

    args = parser.parse_args()

    initial_window = None
    if args.window_init is not None:
        with open(args.window_init, "rb") as f:
            initial_window = list(f.read())

    if args.decompress:
        decoder = LZ77SlidingWindowDecoder(initial_window=initial_window)
        decoder.decode_file(args.input, args.output)
    else:
        # match_finder = MorphingMatchChainBasedMatchFinder(hash_length=1, minimum_match_length=1)
        match_finder = MorphingMatchChainBasedMatchFinder()
        if args.max_chain_length:
            match_finder.max_chain_length = args.max_chain_length
        encoder = LZ77SlidingWindowEncoder(match_finder, initial_window=initial_window)
        encoder.encode_file(args.input, args.output, block_size=BLOCKSIZE)
