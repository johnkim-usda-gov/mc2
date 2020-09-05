	
c*****************************************************************
c	$RCSfile: adjlig.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: adjlig.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 19:56:00  bachelet
cheader
c
c******************************************************************
	subroutine adjlig(oldc,frnew,addc,fractl)

	real      oldc, frnew, addc, fractl

c...Adjust the fraction of lignin in structural C when new material
c...  is added.

c    oldc  = grams C in structural before new material is added
c    frnew = fraction of lignin in new structural material
c    addc  = grams structural C being added

c    fractl comes in as fraction of lignin in structural before new
c           material is added; goes out as fraction of lignin in
c           structural with old and new combined.

c...Local variables
	real      newlig, oldlig


c...oldlig  = grams of lignin in existing residue
	oldlig = fractl * oldc

c...newlig = grams of lignin in new residue
	newlig = frnew * addc

c...Compute lignin fraction in combined residue
	fractl = (oldlig + newlig) / (oldc + addc)

	return
	end
