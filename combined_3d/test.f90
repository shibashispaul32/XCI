character(len=20)::filename,a
do i=1,10
write(a,'(i10)')i
filename = 'tube-' // trim(adjustl(a)) // '.dat'
open(unit=i,file=filename)
write(i,*)i
end do
end
