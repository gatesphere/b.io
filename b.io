// b.io
// small, simple bioinformatics library
// Jacob Peck
// BSD License

// info
bio_const_info := Object clone do(
  complements := list(list("A", "T"), list("C", "G"), 
                      list("G", "C"), list("T", "A")) asMap
  codons := list(
    list("UUU", "F"),
    list("UUC", "F"),
    list("UUA", "L"),
    list("UUG", "L"),
    list("UCU", "S"),
    list("UCC", "S"),
    list("UCA", "S"),
    list("UCG", "S"),
    list("UAU", "Y"),
    list("UAC", "Y"),
    list("UAA", "Stop"),
    list("UAG", "Stop"),
    list("UGU", "C"),
    list("UGC", "C"),
    list("UGA", "Stop"),
    list("UGG", "W"),
    list("CUU", "L"),
    list("CUC", "L"),
    list("CUA", "L"),
    list("CUG", "L"),
    list("CCU", "P"),
    list("CCC", "P"),
    list("CCA", "P"),
    list("CCG", "P"),
    list("CAU", "H"),
    list("CAC", "H"),
    list("CAA", "Q"),
    list("CAG", "Q"),
    list("CGU", "R"),
    list("CGC", "R"),
    list("CGA", "R"),
    list("CGG", "R"),
    list("AUU", "I"),
    list("AUC", "I"),
    list("AUA", "I"),
    list("AUG", "M"),
    list("ACU", "T"),
    list("ACC", "T"),
    list("ACA", "T"),
    list("ACG", "T"),
    list("AAU", "N"),
    list("AAC", "N"),
    list("AAA", "K"),
    list("AAG", "K"),
    list("AGU", "S"),
    list("AGC", "S"),
    list("AGA", "R"),
    list("AGG", "R"),
    list("GUU", "V"),
    list("GUC", "V"),
    list("GUA", "V"),
    list("GUG", "V"),
    list("GCU", "A"),
    list("GCC", "A"),
    list("GCA", "A"),
    list("GCG", "A"),
    list("GAU", "D"),
    list("GAC", "D"),
    list("GAA", "E"),
    list("GAG", "E"),
    list("GGU", "G"),
    list("GGC", "G"),
    list("GGA", "G"),
    list("GGG", "G")
  ) asMap
)

// additions to Sequence (for ease of use)
Sequence countNucleotides := method(
  counts := Map clone
  self foreach(c,
    c = c asCharacter
    if(counts at(c) == nil, counts atPut(c, 0))
    counts atPut(c, counts at(c) + 1)
  )
  counts
)

Sequence transcribeToRNA := method(
  self asMutable replaceSeq("T", "U")
)

Sequence transcribeToDNA := method(
  self asMutable replaceSeq("U", "T")
)

Sequence complement := method(
  str := "" asMutable
  self foreach(c,
    c = c asCharacter
    str appendSeq(bio_const_info complements at(c))
  )
  str
)

Sequence hammingDistance := method(str,
  count := 0
  self foreach(i, c,
    if(str at(i) == c, continue, count = count+1)
  )
  count
)

Sequence triplets := method(
  out := list
  gather := "" asMutable
  index := 0
  self foreach(c,
    c = c asCharacter
    gather appendSeq(c)
    index = index + 1
    if(index == 3,
      out append(gather)
      gather = "" asMutable
      index = 0
    )
  )
  out
)

Sequence toCodonString := method(skipstop,
  str := self triplets
  out := "" asMutable
  str foreach(trip,
    if(bio_const_info codons at(trip) == "Stop" and skipstop,
      continue
      ,
      out appendSeq(bio_const_info codons at(trip))
    )
  )
  out
)

Sequence locations := method(str,
  positions := list
  index := 0
  loop(
    pos := self findSeq(str, index)
    if(pos == nil, break)
    index = pos + 1
    positions append(index)
  )
  positions
)

