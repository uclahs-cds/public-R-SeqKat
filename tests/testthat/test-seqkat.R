expect_error(
	{
	seqkat(5, 3.2, 4, bed.dir = ".")
		},
	"Please supply a path to the reference genome with the ref.dir argument.",
	info = "seqkat Error 1"
	);