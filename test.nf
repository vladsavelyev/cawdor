
line = "8414746.r-man2  vs2870   express- exp16m100   10956   1  16  100gb 24:00 R 19:50"
group = (line =~ /(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+).*/)[1]
println(group[0])
println(group[1])

