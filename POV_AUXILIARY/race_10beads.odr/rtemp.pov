#version 3.6;

#include "colors.inc"
#include "stones1.inc"
#include "stones2.inc"
#include "metals.inc"

global_settings {
	max_trace_level 64
}

#declare zoom=0.105;

camera {
	location <20,-100,80>
	sky z
	right -zoom*x*image_width/image_height
	up zoom*z
	look_at <0,0,0>
}

background{rgb 0}

light_source{<-8,-20,30> color rgb <0.97,0.95,0.95>}
light_source{<25,-12,12> color rgb <0.68,0.70,0.70>}

#declare fv=0;
#declare t_mesh=texture {
	pigment {
		gradient z
		color_map {
			[0 color rgb <0,0,0.6> filter fv]
			[0.15 color rgb <0.2,0.4,1> filter fv]
			[0.5 color rgb <0.6,0.1,0.7> filter fv]
			[0.85 color rgb <1,0.1,0.1> filter fv]
			[1 color rgb <1,0.6,0> filter fv]
		}
		scale <1,1,0.5>
		translate <0,0,-0.35>
	}
	finish {
		ambient 0.34
		diffuse 0.78
		specular 0.24
		phong 0.25
	}
}

#declare f1=finish{F_MetalC}
#declare p0=pigment{rgb <0.97,0.4,0.4>}
#declare p1=pigment{rgb <0.97,0.65,0.32>}
#declare p2=pigment{rgb <0.88,0.88,0.4>}
#declare p3=pigment{rgb <0.4,0.92,0.3>}
#declare p4=pigment{rgb <0.5,0.5,0.97>}
#declare p5=pigment{rgb <0.88,0.42,0.97>}
#declare p6=pigment{rgb <0.97,0.4,0.4>}
#declare p7=pigment{rgb <0.97,0.65,0.32>}
#declare p8=pigment{rgb <0.88,0.88,0.4>}
#declare p9=pigment{rgb <0.4,0.92,0.3>}
#declare p10=pigment{rgb <0.5,0.5,0.97>}
#declare p11=pigment{rgb <0.88,0.42,0.97>}

union {
sphere{<0,0,0>,0.5 texture{T_Stone1 translate <4,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,3.75652,0.870175,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone2 translate <8,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,5.35072,1.77597,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone3 translate <12,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,2.95942,-3.26026,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone4 translate <16,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,1.47391,4.23974,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone5 translate <20,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,4.80725,0.870175,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone6 translate <24,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,2.27101,-5.14432,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone7 translate <28,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,4.87971,-2.71678,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone8 translate <32,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,5.93043,-2.71678,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone9 translate <36,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,3.86522,-0.94142,0.499>}
sphere{<0,0,0>,0.5 texture{T_Stone10 translate <40,0,0>} matrix<1,1.7659e-09,-1.7527e-06,3.5683e-09,1,0.0030433,1.7527e-06,-0.0030433,1,4.87971,-0.94142,0.499>}
    //texture{finish{f1} pigment{p2}}
}

#declare c30=sqrt(0.75);
#declare cr=5.72*1;
#declare th=0.2;
#declare ct=cr+th;

//    cylinder{<1.8,0.00215731,-th>,<1.8,0.00215731,1>,cr+th}
//    cylinder{<1.8,0.00215731,0>,<1.8,0.00215731,1.5>,cr}
prism{
    -th,1,14,
    <cr/c30,0>,<cr/sqrt(3.),cr>,<-cr/sqrt(3.),cr>,<-cr/c30,0>,<-cr/sqrt(3.),-cr>,<cr/sqrt(3.),-cr>,<cr/c30,0>,
    <ct/c30,0>,<ct/sqrt(3.),ct>,<-ct/sqrt(3.),ct>,<-ct/c30,0>,<-ct/sqrt(3.),-ct>,<ct/sqrt(3.),-ct>,<ct/c30,0>
    rotate <90,0,0>
    translate <1.8,0.00215731,0>
    texture{finish{F_MetalA} pigment{rgb <0.2,0.3,0.7>}}
}

#declare ch=cr+0.5*th;

prism{
    -0.5*th,0,7,
    <ch/c30,0>,<ch/sqrt(3.),ch>,<-ch/sqrt(3.),ch>,<-ch/c30,0>,<-ch/sqrt(3.),-ch>,<ch/sqrt(3.),-ch>,<ch/c30,0>
    rotate <90,0,0>
    translate <1.8,0.00215731,0>
    texture{finish{F_MetalB} pigment{rgb <0.5,0.5,0.5>}}
}
