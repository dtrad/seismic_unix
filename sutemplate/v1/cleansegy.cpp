#include "su.h"
#include "segy.h"
segy cleansegy(segy tr)
{
  fprintf(stderr,"cleaning segy..\n");
	tr.tracl=0;
	tr.tracr=0;	/* trace sequence number within reel */
	tr.fldr=0;	/* field record number */
	tr.tracf=0;	/* trace number within field record */
	tr.ep=0;	/* energy source potr.number */
	tr.cdp=0;	/* CDP ensemble number */
	tr.cdpt=0;	/* trace number within CDP ensemble */
	tr.trid=0;	/* trace identification code:*/
	tr.nvs=0;	/* number of vertically summed traces (see vscode
			   in bhed structure) */
	tr.nhs=0;	/* number of horizontally summed traces (see vscode
			   in bhed structure) */
	tr.duse=0;	/* data use:
				1 = production
				2 = test */
	tr.offset=0;	/* distance from source point to receiver
			   group (negative if opposite to direction
			   in which the line was shot) */
	tr.gelev=0;	/* receiver group elevation from sea level
			   (above sea level is positive) */
	tr.selev=0;	/* source elevation from sea level
			   (above sea level is positive) */
	tr.sdepth=0;	/* source depth (positive) */
	tr.gdel=0;	/* datum elevation at receiver group */
	tr.sdel=0;	/* datum elevation at source */
	tr.swdep=0;	/* water depth at source */
	tr.gwdep=0;	/* water depth at receiver group */
	tr.scalel=0;	/* scale factor for previous 7 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */
	tr.scalco=0;	/* scale factor for next 4 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */
	tr. sx=0;	/* X source coordinate */
	tr. sy=0;	/* Y source coordinate */
	tr. gx=0;	/* X group coordinate */
	tr. gy=0;	/* Y group coordinate */
	tr.counit=0;	/* coordinate units code:
				for previous four entries
				1 = length (meters or feet)
				2 = seconds of arc (in this case, the
				X values are longitude and the Y values
				are latitude, a positive value designates
				the number of seconds east of Greenwich
				or north of the equator */
	tr.wevel=0;	/* weathering velocity */
	tr.swevel=0;	/* subweathering velocity */
	tr.sut=0;	/* uphole time at source */
	tr.gut=0;	/* uphole time at receiver group */
	tr.sstat=0;	/* source static correction */
	tr.gstat=0;	/* group static correction */
	tr.tstat=0;	/* total static applied */
	tr.laga=0;	/* lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */

	tr.lagb=0;	/* lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */

	tr.delrt=0;	/* delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */

	tr.muts=0;	/* mute time--start */
	tr.mute=0;	/* mute time--end */
	tr.ns=0;	/* number of samples in this trace */
	tr.dt=0;	/* sample interval=0; in micro-seconds */
	tr.gain=0;	/* gain type of field instruments code:
				1 = fixed
				2 = binary
				3 = floating point
				4 ---- N = optional use */
	tr.igc=0;	/* instrument gain constant */
	tr.igi=0;	/* instrument early or initial gain */
	tr.corr=0;	/* correlated:
				1 = no
				2 = yes */
	tr.sfs=0;	/* sweep frequency at start */
	tr.sfe=0;	/* sweep frequency at end */
	tr.slen=0;	/* sweep length in ms */
	tr.styp=0;	/* sweep type code:
				1 = linear
				2 = cos-squared
				3 = other */
	tr.stas=0;	/* sweep trace length at start in ms */
	tr.stae=0;	/* sweep trace length at end in ms */
	tr.tatyp=0;	/* taper type: 1=linear, 2=cos^2, 3=other */
	tr.afilf=0;	/* alias filter frequency if used */
	tr.afils=0;	/* alias filter slope */
	tr.nofilf=0;	/* notch filter frequency if used */
	tr.nofils=0;	/* notch filter slope */
	tr.lcf=0;	/* low cut frequency if used */
	tr.hcf=0;	/* high cut frequncy if used */
	tr.lcs=0;	/* low cut slope */
	tr.hcs=0;	/* high cut slope */
	tr.year=0;	/* year data recorded */
	tr.day=0;	/* day of year */
	tr.hour=0;	/* hour of day (24 hour clock) */
	tr.minute=0;	/* minute of hour */
	tr.sec=0;	/* second of minute */
	tr.timbas=0;	/* time basis code:
				1 = local
				2 = GMT
				3 = other */
	tr.trwf=0;	/* trace weighting factor, defined as 1/2^N
			   volts for the least sigificant bit */
	tr.grnors=0;	/* geophone group number of roll switch
			   position one */
	tr.grnofr=0;	/* geophone group number of trace one within
			   original field record */
	tr.grnlof=0;	/* geophone group number of last trace within
			   original field record */
	tr.gaps=0;	/* gap size (total number of groups dropped) */
	tr.otrav=0;	/* overtravel taper code:
				1 = down (or behind)
				2 = up (or ahead) */
	/* local assignments */
	tr.d1=0;	/* sample spacing for non-seismic data */
	tr.f1=0;	/* first sample location for non-seismic data */
	tr.d2=0;	/* sample spacing between traces */
	tr.f2=0;	/* first trace location */
	tr.ungpow=0;	/* negative of power used for dynamic
			   range compression */
	tr.unscale=0;	/* reciprocal of scaling factor to normalize
			   range */
	tr.ntr=0; 	/* number of traces */
	tr.mark=0;	/* mark selected traces */

        tr.shortpad=0; /* alignment padding */
	return(tr);
}




