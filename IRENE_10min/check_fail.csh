#! /bin/csh

### Transects ==> j loop
set j=0
set s=0
set f=0
while($j<4158)
##echo ==> print, must use $ before variable
  if ( -d T$j) then 
    cd T$j/
    set File = XBeach_CR*
    grep -q 'Successfully completed.' $File
    if ( $status == 0 ) then
      #echo "Success: T$j"
      @ s++
    else
      echo "Fail: T$j"
      @ f++
    endif
    cd ..
  else
    #echo "no T$j"
  endif
  @ j++
end

echo "Successful: $s"
echo "Failed: $f"
