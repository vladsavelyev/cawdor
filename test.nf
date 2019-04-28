params.str = 'Hello world!'

process splitLetters {

  output:
  file 'chunk_*' into letters mode flatten

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

letters = letters.view {
  "  Letter: ${it.baseName}"
}