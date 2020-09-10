      
c*****************************************************************
c	$RCSfile: eq_test.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: eq_test.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:20:24  bachelet
cheader
c
c******************************************************************
	subroutine eq_test(eq_flag, min_years, final_growth)

c-----C. Daly - 2/23/96.  For standard Century runs, run model only until
c-----slow pool organic matter (som2c) achieves relative equilibrium.

c-----C. Daly 5-6-96.  Split off of cen_step.

c-----C. Daly 6-11-96.  Now bringing in eqyears, read in place of
c-----isotope labeling type in sch file header (detiv).  Eqyears gives 
c-----disturbance cycle over which eq code should test for som2c eq.
c-----------------------------------------------------------------------
      include 'zztim.inc'
      include 'schvar.inc'
      include 'timvar.inc'
      include 't0par.inc'
      include 'plot1.inc'
      include 'mapssvar.inc'

      logical eq_flag
      logical eq_start
      save eq_start
      real som2c_last, threshold, final_growth, athresh
      save som2c_last
      integer min_years
      integer counter
      save counter

c-----Set EQ thresholds for som2c change in fraction per year.
      data threshold /0.0001/, athresh /0.01/

c-----Set fraction of last cycle to simulate.
c-----For one-step init, use 0.0
c-----For normal runs, use 0.0 or whatever.
c-----For last-year output, use 0.67 or so.
c*****Now passed in as a mapss switch (-g)
cx      data final_growth /0.0/
c---------------------------------------------------------------------



c-----Initialize stuff first time around.
      if( nint(time) .eq. 1 .and. month .eq. 1) then
c-------Year counter for getting part way into last cycle.
        counter = 0
c-------EQ flag to "not in equilibrium".
        eq_flag = .false.
        eq_start = .false.
c-------Check eqyears; set to default if zero.
        if(eqyears .eq. 0) eqyears = rptyrs
        return
      endif

c-----If final growth counter is ticking, check for EQ time.
      if (eq_start) then 
        if (counter .ge. float(eqyears) * final_growth * 12
     .      .and. month .eq. 12) then
          eq_flag = .true.
          if (diags) write(*,*) 'Final time = ',time
          return
        else
          counter = counter + 1
          return
        endif
      endif

c-----Never check EQ in the middle of the schedule cycle or other than in
c-----December.
      if (month .ne. 12) return
      if (mod (nint(time - 1.), eqyears) .ne. 0) return
        
c-----Initialize som2c at end of very first cycle. 
      if ( nint(time - 1.) .eq. eqyears) then
c-------Last time's som2c value
        som2c_last = som2c(1)
        return
      endif

c-----Run for at least min_years or for one cycle, whichever is
c-----greatest, to get things going properly.
      if (time .lt. max(min_years, eqyears) ) return

      if (diags) then
         write(*,*) 'Slow pool soil organic matter: check time = ',
     .           time, ' som2c = ',som2c
      endif

c-----If time is not up, and change in som2c is still large, continue on.
      if ((tend - eqyears - time) .lt. (dt * 0.5)) then
        if (diags) write(*,*) 'Ding!  Century says times up, sodbuster!'
c        eq_flag = .true.
c        return
        eq_start = .true.
      else
c-----If som2c pool is stable, initiate counter and proceed part
c-----way through next cycle before exiting.  This will sample veg
c-----in a "middle" phase of development for lai maps, etc. 

c.....Thinking we might get a problem if som2c becomes zero before
c.....equilibrium is reached so adding another condition here
c.....D. Bachelet 10-8-97

	if(som2c(1) .eq. 0.0) then
		eq_start = .true.
		eq_flag = .true.
		if (diags) write(*,*) 'Dong! som2c = zero'
		return
	endif 
      change = abs(som2c(1) - som2c_last)
c      write(*,*) 'time, som2c(1), change', time, som2c(1), change
        if ((abs(som2c(1) - som2c_last) .gt. threshold*eqyears*som2c(1))
     .    .and. (abs(som2c(1) - som2c_last) .gt. athresh)) then
          som2c_last = som2c(1)
          return
        else
          eq_start = .true.
c---------Check and see if we should stop now.  If final_growth parameter
c---------is less than 1 year's worth of time into the next disturbance
c---------cycle, end it.
          if (final_growth*eqyears .lt. 1.0) then
            eq_flag = .true.
            if (diags) write(*,*) 'som2c equilibrium at ',time,
     .               ' - No Post-EQ cycle needed.'
          else
            if (diags) write(*,*) 'som2c equilibrium at ',time,
     .               ' - counter is running...'
          endif
        endif
      endif

      return
      end
