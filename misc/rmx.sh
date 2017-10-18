rmx1() {
if [ -e $1 ]; then

mkdir ../temporaryfolder
 cp -r $1 ../temporaryfolder
 rm *
 cp -r ../temporaryfolder/* .
 rm -r ../temporaryfolder
fi
}

rmx2() {
if [ -e $1 ] && [ -e $2 ];  then

mkdir ../temporaryfolder
 cp -r $1 ../temporaryfolder
 cp -r $2 ../temporaryfolder
 rm  *
 cp -r ../temporaryfolder/* .
 rm -r ../temporaryfolder
fi
}

rmx3() {
if [ -e $1 ] && [ -e $2 ] && [ -e $3 ]; then

mkdir ../temporaryfolder
 cp -r $1 ../temporaryfolder
 cp -r $2 ../temporaryfolder
 cp -r $3 ../temporaryfolder
 rm  *
 cp -r ../temporaryfolder/* .
 rm -r ../temporaryfolder
fi
}

if [ "$#" -eq 3 ]; then
 rmx3 $1 $2 $3
elif [ "$#" -eq 2 ]; then
 rmx2 $1 $2
elif [ "$#" -eq 1 ]; then
 rmx1 $1
else
 echo "Bad # of arguments"
fi
