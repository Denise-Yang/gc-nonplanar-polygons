// The MIT License
// Copyright © 2021 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Exact distance to a sphere cut by a plane. Beware doing the max() of
// a sphere and a plane won't produce an exact Euclidean distance.
// Based on sdCutDisk(): https://www.shadertoy.com/view/ftVXRc
//
// It is a useful primitive when combined with rounding/inflating, which
// cannot be done with the non-Euclidean max() approach, since you can do
// things like mushroom heads.

// List of other 3D SDFs: https://www.shadertoy.com/playlist/43cXRl
//
// and https://iquilezles.org/articles/distfunctions



const vec3 sphereCenter = vec3(0.f, 0.f, 0.f);
const float outerRadius = 1.25f;

//const float shift = -.25*(-1.f - 4.f * outerRadius*outerRadius);
//const float shift =  outerRadius*outerRadius + outerRadius;




//when changing the spherical harmonic function, make sure to update, sphericalHarmonic(), map(), and gradShade()

/** Start Helper Functions**/
float getDis(vec3 pt1, vec3 pt2){
  return sqrt(pow(pt1[0]-pt2[0],2.f)  + pow(pt2[1]-pt2[1],2.f) + pow(pt1[2]-pt2[2],2.f));

}

float getDis2D(vec2 pt1, vec2 pt2){
  return sqrt(pow(pt1[0]-pt2[0],2.f)  + pow(pt2[1]-pt2[1],2.f));

}


//minimize function of form x^2 - z^2 + y + c inside ball of radius p
float findShift(float p){
    return -.25*(-1.f - 4.f * p*p);
}





//For harnack's checks if we're close to our desired level set
bool closeToLevelset(float ang, float levelset){
    float epsilon = .001;
    return ang < levelset + epsilon  && ang > levelset - epsilon;
}

float getMaxStep(float fx, float R, float lo_bound, float up_bound){
    float v = up_bound/fx;
    float w = lo_bound/fx;
    float lo_r = (-(2.*w+1.) + sqrt(8.*v + 1.))*R/(2.*w);
    float up_r = (2.*v+1. - sqrt(8.*v + 1.))*R/(2.*v);
    if (up_r < lo_r) return up_r;
    return lo_r;
}

float getMaxStep2D(float fx, float R, float lo_bound, float up_bound, float v0, float vz){
    float lo_bound_b =fx+v0+R*vz;
    float up_bound_b =fx+v0-R*vz;
    float a = vz;
    float c = -v0*R + fx*R;
    float lo_r = (-lo_bound_b + sqrt(lo_bound_b*lo_bound_b-4.*a*c))/(2.*a);
    float up_r = (-up_bound_b - sqrt(up_bound_b*up_bound_b-4.*a*c))/(2.*a);
    if (up_r < lo_r) return up_r;
    return lo_r;
}


float sphericalHarmonic(vec3 pt, float shift){
  float x = pt[0];
  float y = pt[1];
  float z = pt[2];
  return pow(x,2.f) - pow(z,2.f) - y + shift;
}

float sphericalHarmonic2D(vec2 pt, float shift){
  float x = pt[0];
  float z = pt[1];
  return pow(x,2.f) - pow(z,2.f) + shift;
}


bool outsideCircle( in vec2 pos )
    
{
    float x = pos[0];
    float y = pos[1];
    return (x*x + (y-sphereCenter[1])*(y-sphereCenter[1]) - 1.f > 0.f);
}


bool outsideSphere( in vec3 pos )
    
{
    float x = pos[0];
    float y = pos[1];
    float z = pos[2];
    return (x*x + (y-sphereCenter[1])*(y-sphereCenter[1]) + z*z - 1.f > 0.f);
}

float map( in vec3 pos )
    
{
    float x = pos[0];
    float y = pos[1];
    float z = pos[2];
    if (outsideSphere(pos)){
        return 1.f;
    }
    return x*x - z*z + y + 1.f;
}



//something buggy about this function or like we're not getting accurate harnack results because of it
bool intersectSphere(vec3 ro, vec3 rd, inout float t0, inout float t1){
  float t = -1.0f;
  vec3 L = ro - sphereCenter;
  float a = dot(rd,rd);
  float b = 2.0f * dot(rd,L);
  float c = dot(L,L)-1.f;
  float discr = b * b - 4.f * a * c;
  if (discr < 0.f) return false;
  else if (discr == 0.f){
    t0 = t1 = - 0.5 * b / a;
    return true;
  }
  else{
    float q = (b > 0.f) ?
            -0.5 * (b + sqrt(discr)) :
            -0.5 * (b - sqrt(discr));
        t0 = q / a;
        t1 = c / q;
        if (t1 < t0) {
            float temp = t0;
            t0 = t1;
            t1 = temp;
        }
        //might not want to return t1 instead since then the ray would be hittingthe edge of the sphere
       
        return true;
  }
}

/** End Helper Functions**/



// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
    vec2 e = vec2(1.0,-1.0)*0.5773;
    const float eps = 0.0005;
    return normalize( e.xyy*map( pos + e.xyy*eps ) + 
					  e.yyx*map( pos + e.yyx*eps ) + 
					  e.yxy*map( pos + e.yxy*eps ) + 
					  e.xxx*map( pos + e.xxx*eps ) );
}

vec3 gradShade( in vec3 pos )
{
    float x = pos[0];
    float y = pos[1];
    float z = pos[2];
    float mag = sqrt(4.f*x*x + 1.f + 4.f*z*z);
    return vec3(x/mag, y/mag, z/mag);
}





#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 3
#endif

bool raymarch(vec3 ro, vec3 rd, float tmax, inout vec3 finalPos){
    float t = 0.0;
    for( int i=0; i<350; i++ )
        {
            vec3 pos = ro + t*rd;
            float h = map(pos);
            if( h<0.0001 || t>tmax ) break;
            t += 0.05;
        }
    finalPos = ro + t*rd;
    return (t < tmax);
        
}

bool harnack(vec3 ro, vec3 rd,  inout vec3 pos, inout bool maxSteps, float time){
    int i = 0;
    
    float t0;
    float t1;
    float tmax;
    float t = 0.f;
    float c = 0.f;
    float domainRadius = 1.f;
    //float shift = findShift(domainRadius);
    float shift = 1.25*1.25;
    float levelset = 0.f+ shift;
    maxSteps = false;
    
    vec2 ro2 = ro.xz;
    vec2 rd2 = rd.xz;
    
    
    
    
        
        
    bool hitSphere = intersectSphere(ro, rd, t0, t1);
    if (!hitSphere){
        return false;
    }
    
    if (t0 < 0.f && t1 < 0.f) return false;
    if (t0  > 0.f) t = t0;
    tmax = t1;
    
    
    pos = ro+t*rd;
    vec2 pos2 = pos.xz;
    {    
        float loBound = -.43;
        float hiBound = levelset+pos.y;
        if (sphericalHarmonic2D(pos2, shift) > levelset){
            loBound = levelset+pos.y;
            hiBound = 10.f;
        }
        while (t < tmax){
            if (i > 10000){
                maxSteps = true;
                return false;
            }
            float val = sphericalHarmonic2D(pos2, shift); 
            if (closeToLevelset(val, shift+pos.y)){
                return true;
            }
            
            //float R = outerRadius - getDis(pos,sphereCenter);
          
            float R = outerRadius - getDis2D(pos2,sphereCenter.xz);
            //float r = getMaxStep(val, R, loBound, hiBound);
            float r = getMaxStep2D(val, R, loBound, hiBound, pos[1],rd[1]);
            t += r;
            pos = ro + t*rd;
            pos2 = pos.xz;
            i++;
        }
        
     }
     return false;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //toggle harnacks
    bool doHarnacks = true;
    bool maxSteps = false;
    
    // camera movement	
	float an = 0.8*iTime;
	vec3 ro = vec3( 3.0*cos(an), -.5f, 3.0*sin(an) );
	//vec3 ro = vec3( 2., 0.0, 2. );
    vec3 ta = vec3( 0.0, 0.0, 0.0 );
    // camera matrix
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    float t;
    
    float f1,f2;
    
    
    vec3 tot = vec3(0.0);
    
    #if AA>1
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (2.0*(fragCoord+o)-iResolution.xy)/iResolution.y;
        #else    
        vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
        #endif

	    // create view ray
        vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

        // raymarch
        const float tmax = 5.0;
        vec3 pos;
        bool didHit = false;
        if (doHarnacks){
            didHit = harnack(ro,rd, pos, maxSteps, iTime);
        } else {
            didHit = raymarch(ro,rd,tmax, pos);
        }
    
        // shading/lighting	
        vec3 col = vec3(0.0);
        if (maxSteps){
            col =  vec3(1.f,0.0,0.0);
        }
        else if( didHit )
        {
            vec3 nor;
            if (doHarnacks){
               nor = gradShade(pos);
            } else{ 
               nor = calcNormal(pos);
            }
            float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0));
            col = vec3(0.2,0.3,0.4)*amb + vec3(0.8,0.7,0.5)*dif;
        }else if (intersectSphere(ro,rd, f1,f2)){
        if (f1  > 0.f) t = f1;
            else{
                t = f2;
            }
            
            vec3 nor = gradShade(ro+rd*t);
            float amb = 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0));
            
            col = vec3(0.5,0.1,0.1)*amb;
        }
        

        // gamma        
        col = sqrt( col );
	    tot += col;
    #if AA>1
    }
    tot /= float(AA*AA);
    #endif

	fragColor = vec4( tot, 1.0 );
}
