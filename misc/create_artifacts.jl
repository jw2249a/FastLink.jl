using SHA, Inflate, Tar
fileA = "test/dfA.tar.gz"
fileB = "test/dfB.tar.gz"

println("sha256 = ", bytes2hex(open(sha256, fileA)))
println("git-tree-sha1 = ", Tar.tree_hash(IOBuffer(inflate_gzip(fileA))))

println("sha256 = ", bytes2hex(open(sha256, fileB)))
println("git-tree-sha1 = ", Tar.tree_hash(IOBuffer(inflate_gzip(fileB))))
 
