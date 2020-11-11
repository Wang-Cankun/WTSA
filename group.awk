$1 != f1 { if (NR>1) print f2; f1=f2=s="" }
{ f1=$1; f2=f2 s $2; s="," }
END { print f2 }
