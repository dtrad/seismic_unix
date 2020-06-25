#! /bin/sh

# setenv CWPROOT /home/tkrishna/SU
# setenv PATH $CWPROOT/bin":$PATH"

#===================================================
# SCRIPT TO GENERATE MODELS AND SHOOT A SINGLE
# SHOT TO GENERATE BOTH PP AND PS SEISMIC DATA
# THIS WILL ALSO GENERATE A 2C VSP DATA
#
# THE MODEL LAYERS IS DEFINED IN LAYERS SHOWN IN 
# THE FILE model001.unif
# THE FORMAT IS:
# EACH LAYER IS SEPARATED BY 1 -99999
# THE FIRST COLUMN IS THE X AND THE SECOND IS THE Z 
#
#    0         0
# 5000         0
#    1    -99999
#    0       100
# 5000       100
#    1    -99999
#    0       350
# 5000.      450
#    1    -99999
#    0       600
# 5000.      600
#    1.   -99999
#    0       850
# 5000       850
#    1    -99999
#    0      1000
# 5000      1000
#    1    -99999
#
#===================================================

# Select what parts of this script you would like to run

makemodel=true
plotmodel=false
runfd=false
plotseismic=false

# Define parameters for the modeling and output
modelfile=model001.unif         # input model file for unif2aniso

ninf=4                          # number of interfaces (surface counts)

x0=0,0,0,0,0                    # x-position(s) for  vp00,vs00,rho00, etc.
z0=0,5,10,15,20			# z-position(s) for  vp00,vs00,rho00, etc.
nz=400                          # size of z (depth) dimension  of model
nx=1000                         # size of x (horizontal) dimension of model
dz=2.5                          # increment in z direction (m)
dx=5                            # increment in x direction (m)

vp00=1500,1800,2200,2600,3000   # P-wavespeed(s) at (z0,x0)
vs00=0,500,950,1300,1800        # S-wavespeed(s) at (z0,x0)
kz=0.0,1.2,0.5,0.0,0.0          # z-derivative of model (dv*/dz)
rho00=1000,1200,1700,2100,2400  # density(s) at (z0,x0)
eps00=0.0,0.0,0.0,0.0,0.0       # thompson's epsilon
del00=0.0,0.0,0.0,0.0,0.0       # thompson's delta
q00=0.0,0.0,0.0,0.0,0.0         # Q factor

method=mono			# boundary interpolation scheme for unif2aniso
                                #     linear=interpolation of interface 
                                #     =mono for monotonic cubic,=akima for Akima's cubic interpolation
                                #     =spline for cubic spline interpolation
                                
dt=0.0005                       # time sampling interval (s)
lt=2.5	                        # latest time modeled (s)
fx=0                            # first x value
verbose=1                       # =1 chatty, =0 silent
snfile="snaps.su"               # output file for snapshots
rhofile="rho_file"              # input file of densities
snaptime=0.02,0.05,.1,.15,.2,.25.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95	# times of snapshots

bc=2,1000,1000,1000             # boundary conditions, Top,left,bottom,right 
                                # (=0 none, =1 symmetry,=2 free surface (top only)
                                #  >2 absorbing (value= width of absorbing layer)

qsw=0                           # =1 put in attenuation
asw=0                           # =1 anisotropy

favg=25.0                       # average frequency
ts=.05                          # source duration
wtype=dg                        # waveform type, dg: Gaussian derivative , ga: Gaussia, ri: Ricker , sp: spike, sp2: double spike

# shot & receiver position
hsz=100		                # z-position of horizontal line of geophones
vsx=500		                # x-position of vertical line of geophones (VSP)
sx=1475                         # x-position of sources
sz=10                           # z-position of sources

# output directory
outdir="/data/cgyv90b_vd1/tkrishna/ALBACORA/SYNTH002"

# plot parameters
cmap=hsv2  #which colour map to use (R-G-B)

if($makemodel)
then

# build stiffness and density files
unif2aniso < $modelfile ninf=$ninf x0=$x0 z0=$z0 nz=$nz nx=$nx \
dx=$dx dz=$dz vp00=$vp00 vs00=$vs00 rho00=$rho00 eps00=$eps00 delta00=$del00 method=$method

# this is just to make some QC items to have a look
# at the velocities and anisotropic models
unif2 < $modelfile ninf=$ninf nz=$nz nx=$nx \
dz=$dz dx=$dx v00=$vp00 dvdz=$kz method=$method >vp_file

unif2 < $modelfile ninf=$ninf nz=$nz nx=$nx \
dz=$dz dx=$dx v00=$vs00 dvdz=$kz method=$method >vs_file

unif2 < $modelfile ninf=$ninf nz=$nz nx=$nx \
dz=$dz dx=$dx v00=$eps00 method=$method >eps_file

unif2 < $modelfile ninf=$ninf nz=$nz nx=$nx \
dz=$dz dx=$dx v00=$del00  method=$method >del_file

xbox=10
ybox=10
nxplot=`bc -l <<-END
        scale=0
        $nx  * 2
END`
nzplot=`bc -l <<-END
        scale=0
        $nz  * 2
END`
fi

# Plot the models to screen
if($plotmodel)
then
for i in vp vs  rho # del eps c11 c13 c15 c33 c35 c55
do
	echo $xbox $ybox

	ximage <  ${i}_file n1=$nz d1=$dz n2=$nx d2=$dx legend=1 \
		title=" ${i} parameter file" cmap=$cmap	&

	xbox=`expr $xbox + 110 `
	ybox=`expr $ybox + 5 `

done

fi


# transpose stiffness and density for modeling
for i in c11 c13 c15 c33 c35 c55 rho
do
	echo "Transposing Stiffness Models"

        mv ${i}_file tmp.file
        transp n1=$nz < tmp.file > ${i}_file
done
rm tmp.file

# Run the finite difference modeling
if($runfd)
then

for sx in {2000..2050..25}
do
suea2df dt=$dt lt=$lt nz=$nz fx=$fx nx=$nx dx=$dx dz=$dz verbose=1 \
rhofile=$rhofile hsz=$hsz vsx=$vsx \
bc=$bc qsw=$qsw asw=$asw sx=$sx sz=$sz favg=$favg ts=$ts wtype=$wtype \
vrsfile=$outdir/vsp_${sx} hsfile=$outdir/hs_${sx}
done
fi


# Plot the seismic to screen
if($plotseismic)
then
n2=`bc -l <<-END
      scale=1
       $nx * 2
END`

echo $n2


#suxmovie < snaps.su n1=$nz n2=$n2 clip=1e-12 loop=1 title="snapshots horizontal vertical " width=$nxplot height=$nzplot  sleep=200000 &

for sht in {2000..2050..25}
do
# shot gathers from a horizontal line of geophones
suximage <  $outdir/hs_${sht}.su xbox=0 ybox=400 wbox=$nxplot hbox=$nzplot  perc=99 title="Shot gathers  Vertical Horizontal " &

# this is to plot the VSP
#suximage <  $outdir/vsp_${sht}.su xbox=400 ybox=400  wbox=$nxplot hbox=$nzplot  perc=99 title=" VSP  Vertical Horizontal " &
done

fi
exit 0
