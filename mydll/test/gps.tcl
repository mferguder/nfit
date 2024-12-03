namespace import blt::*


proc readfile {fid} {
    while { [gets $fid line] > 3 } {
	set x [lindex $line 0]
	set y [lindex $line 1]
	set rx [lindex $line 2]
	set ry [lindex $line 3]
	set x [expr -$x]
	push $x $y $rx $ry
    }
}

set thresh 100

proc transfile {fid fout a b c d} {
    global thresh
    set iter 0
    while { [gets $fid line] > 3 } {
	set x [lindex $line 0]
	set y [lindex $line 1]
	set rx [lindex $line 2]
	set ry [lindex $line 3]
#	set x [expr -$x]
	set nx [expr -$x*$a-$y*$b+$c]
	set ny [expr $y*$a-$x*$b+$d]
	set dist [expr ($nx-$rx)*($nx-$rx)+($ny-$ry)*($ny-$ry)]
	if { $dist < $thresh } {
	    puts $fout "$x $y $rx $ry"
	   # incr iter
	}
    }
}

proc transfile2 {fid fout a b c d} {
    global thresh
    set sum 0
    while { [gets $fid line] > 3 } {
	set x [lindex $line 0]
	set y [lindex $line 1]
	set rx [lindex $line 2]
	set ry [lindex $line 3]
#	set x [expr -$x]
	set nx [expr -$x*$a-$y*$b+$c]
	set ny [expr $y*$a-$x*$b+$d]
	set dist [expr ($nx-$rx)*($nx-$rx)+($ny-$ry)*($ny-$ry)]
	set sum [expr $sum+$dist]
#	if { $dist < $thresh } {
	    puts $fout "$nx $ny $rx $ry"
#	}
    }
    puts $sum
}

set para "-0.568788 -0.822484 79.8158 -5.66117"


proc step1 {} {
    global para
    set fid [open /home/yufeng/seif3/nf2/path.1 r]
    set fout [open /home/yufeng/seif3/nf2/path.2 w]
    eval transfile $fid $fout $para
    close $fid
    close $fout
}

proc step2 {} {
    set fid [open /home/yufeng/seif3/nf2/path.2 r]
    matcher_init 10000
    readfile $fid
    close $fid
    return [match]
}

proc step3 {} {
    global para
    set fid [open /home/yufeng/seif3/nf2/path.2 r]
    set fout [open /home/yufeng/seif3/nf2/path.4 w]
    eval transfile2 $fid $fout $para
    close $fid
    close $fout
}

set opt 0

proc track {fid fout a b c d} {
    global thresh opt
    set iter 0
    while { [gets $fid line] > 3 } {
	set x [lindex $line 0]
	set y [lindex $line 1]
	set rx [lindex $line 2]
	set ry [lindex $line 3]
#	set x [expr -$x]
	set nx [expr -$x*$a-$y*$b+$c]
	set ny [expr $y*$a-$x*$b+$d]
	set dist [expr ($nx-$rx)*($nx-$rx)+($ny-$ry)*($ny-$ry)]

	.g marker configure g -coords "$nx $ny "
	.g marker configure r -coords "$rx $ry "
	update
	vwait opt
	if { $opt == 1 } {
	    if { $dist > 50 } {gets stdin iter}
	    puts $fout "$x $y $rx $ry"
	    puts "$x $y $rx $ry"
	}
    }
}


proc step4 {} {
    global para
    set fid [open /home/yufeng/seif3/nf2/path.1 r]
    set fout [open /home/yufeng/seif3/nf2/path.4 w]
    eval track $fid $fout $para
    close $fid
    close $fout
}

wm geometry . =400x400+0+0
graph .g -width 1024 -height 1024
.g axis configure y -subdivisions 5 -max 200 -min -200 -descending 1
.g axis configure x -subdivisions 5 -max 200 -min -200
pack .g

.g marker create text -coords { 0 0 } -name g -text g -outline red
.g marker create text -coords { 0 0 } -name r -text r -outline green

set zoomMod "Control-"

foreach graph  { .g }  {
    Blt_ZoomStack $graph
    Blt_Crosshairs $graph
}

bind .g <2> {set opt 0; puts $opt}
bind .g <B1-Motion> {set opt 1; puts $opt}
bind .g <1> {set opt 1; puts $opt}
bind .g <B2-Motion> {set opt 0; puts $opt}
