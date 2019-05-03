cd ~/git/trimal/
rm -rf /tmp/new /tmp/old
bin/new_trimal -out /tmp/new -htmlout /tmp/htmlNew.svg $@
../Strimal/source/trimal -out /tmp/old -htmlout /tmp/htmlOld $@
echo ""
# meld /tmp/new /tmp/old &
firefox /tmp/htmlNew.svg /tmp/htmlOld
