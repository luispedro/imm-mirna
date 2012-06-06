import sys

if __name__ == '__main__':
        if len(sys.argv) < 2:
                print "Syntax: python remraw.py <sequence1> <sequence2> ..."
                exit(0)

        inf = sys.stdin
	seqs = sys.argv[1:]

	def markedForRemoval(seq):
		for s in seqs:
			if s == seq:
				return True
		return False

	while True:
		entry = [ inf.readline(), inf.readline(), inf.readline(), inf.readline() ]

		if not entry[0]:
			break

		if not markedForRemoval(entry[1][:-1]):	# remove newline character
			print entry[0], entry[1], entry[2], entry[3],

