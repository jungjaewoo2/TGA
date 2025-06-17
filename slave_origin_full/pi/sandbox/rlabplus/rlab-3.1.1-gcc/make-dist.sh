#
# $1 CVS tag for distribution
# $2 Directory name for distrubution
#
# Must execute from /u1/ian/rlab/rlab2

echo "Use tag:" $1
echo "Create directory:" $2
cvs tag $1
cd ../
cvs export -r $1 -d $2 rlab/rlab2
tar cvf $2.tar $2
gzip -v $2.tar
