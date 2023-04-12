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


const vec3 sphereCenter = vec3(0.f,0.f,0.f);
const float outerRadius = 1.25;

const float M_PI = 3.1415;
const float numPoints = 13.f;
const int len = int(numPoints);
const float numTri = 12.f;
const int len2 = int(numTri);
//vec3 pts[6] = vec3[6]( vec3(1.,0.,1.), vec3(0.,1.,0.), vec3(-1.,0.,.5),  vec3(-.5,.8,.6),vec3(0.,-1.,0.),vec3(0.2,.0,0.) );
//vec3 tri[5] = vec3[5]( vec3(0,1,5), vec3(1,2,5), vec3(2,3,5), vec3(3,4,5), vec3(4,0,5));
/*
vec3 pts[9] = vec3[9]( vec3(.25,0.5,0.), 
                        vec3(0.16667,.16667,0.3333), 
                        vec3(.5,0.,.0),  
                        vec3(.83333,.33333,.0),
                        vec3(0.75,.25,0.5),
                        vec3(0.3333,.3333,0.6667),
                        vec3(0.5,.75,0.5),
                        vec3(0.6667,.6667,0.),
                        vec3(0.4,0.4,0.2)
                        );                      

vec3 pts[9] = vec3[9]( vec3(1.,2.,0.), 
                        vec3(0.6668,.6668,1.333), 
                        vec3(2.,0.,.0),  
                        vec3(3.333,1.33333,.0),
                        vec3(3.,1.,2.),
                        vec3(1.3333,1.3333,2.6667),
                        vec3(2.,3.,2.),
                        vec3(2.6668,2.6668,0.),
                        vec3(0.,0.,0.)
                        );
vec3 tri[8] = vec3[8]( vec3(0,1,8), vec3(1,2,8), vec3(2,3,8), vec3(3,4,8), vec3(4,5,8),vec3(5,6,8),vec3(6,7,8),vec3(7,0,8));

                        */
                        
vec3 pts[13] = vec3[13]( vec3(.833,.1667,0.333), 
                        vec3(0.875,.25,0.0), 
                        vec3(.5,0.,.0),  
                        vec3(.125,.125,.25),
                        vec3(0.1667,.333,0.0),
                        vec3(0.75,.5,0.),
                        vec3(0.6667,.3333,0.6667),
                        vec3(0.375,.375,0.25),
                        vec3(0.5,0.75,0.5),
                        vec3(0.625,.25,0.0),
                        vec3(0.333,0.6667,0.),
                        vec3(0.25,0.25,0.5),
                        vec3(0.0,0.0,0.)
                        );
//vec3 tri[12] = vec3[12]( vec3(0,1,12), vec3(1,2,12), vec3(2,3,12), vec3(3,4,8), vec3(4,5,8),vec3(5,6,8),vec3(6,7,8),vec3(7,0,8), vec3(4,5,8),vec3(5,6,8),vec3(6,7,8),vec3(7,0,8));



//vec3 tri[2] = vec3[2]( vec3(0,1,2), vec3(0,2,3));

//const float shift = -.25*(-1.f - 4.f * outerRadius*outerRadius);
//const float shift =  outerRadius*outerRadius + outerRadius;

vec2 makeComplex(float angle){
    return vec2(cos(angle/2.f), sin(angle/2.f));
}


/** Start Helper Functions**/
float getDis(vec3 pt1, vec3 pt2){
  return sqrt(dot(pt1-pt2, pt1-pt2));

}

float closestPointOnLine(vec3 pt){
    float minD = 100000.f;
    for (int i = 0; i < len; i++){
        vec3 p1 = pts[i];
        vec3 p2 = pts[(i+1)%len];
        vec3 m = p2-p1;
        vec3 v = pt-p1;
        float t = clamp(dot(m,v)/dot(m,m),0.0,1.0); //dot = |a|*|b|cos(theta) * n, isolating |a|sin(theta)
        float d = getDis(v, t * m);

        //float d = length(cross(m,v))/length(m); //cross = |a|*|b|sin(theta) * n, isolating |a|sin(theta)
        
        //get distance of endpoints to get closest dis on line segment
        /*
        float d1 = getDis(pt, p1);
        float d2 = getDis(pt, p2);
        float lineDis = d1 < d2 ? d1 : d2;
        if (lineDis <= d){
            d = lineDis;
        }*/
        if (d < minD){
         minD = d;
        }        
     }
     return minD;

}

float getClosestDis(vec3 pt){
    float minDis = 100000.f;
    for (int i =0;i<int(numPoints);i++){
        float dis = getDis(pts[i],pt);
        if (dis < minDis){
            minDis = dis;
        }
    }
    return minDis;
    
}

float det3( in vec3 a, in vec3 b, in vec3 c) {
    return dot(a, cross(b, c));
}

float inPlane(vec3 pt){
    vec3 n = cross(pts[1]-pts[0], pts[2]-pts[0]);
    return dot(n, pt-pts[0]);
}

bool rayintersect( in vec3 ro, in vec3 rd, out float t ) {
    vec3 a = pts[0] - ro;
    vec3 b = pts[1] - ro;
    vec3 c = pts[2] - ro;
    vec3 r = rd;
    
    float vol = det3(a,b,c);
    float la = det3(r, b, c) / vol;
    float lb = det3(a, r, c) / vol;
    float lc = det3(a, b, r) / vol;
    t = 1.f / (la + lb + lc);
    
    return t >= 0.f && la >= 0.f && lb >= 0.f && lc >= 0.f;
}


bool intersectTriangle(vec3 ro, vec3 rd, inout float intersect){
    vec3 s1 = pts[0] - pts[1];
    vec3 s2 = pts[2] - pts[1];
    vec3 N = cross(s1,s2);
    float NdotDir = dot(N,rd);
    
    if (NdotDir < .01 || NdotDir > -.01){
     return false;
    }
    
    float d = dot(N, pts[0]);
    //float t = -(N.dotProduct( + d) / NdotRayDirection;
    return false;
    

}


float getMaxStep(float fx, float R, float lo_bound, float up_bound){
    float v = up_bound/fx;
    float w = lo_bound/fx;
    float lo_r = (-(2.*w+1.) + sqrt(8.*w + 1.))*R/(2.*w);
    float up_r = (2.*v+1. - sqrt(8.*v + 1.))*R/(2.*v);
    if (up_r < lo_r) return up_r;
    return lo_r;
}


float triangleSA(vec3 p0, vec3 p1, vec3 p2, vec3 p){
    vec3 a = p0-p;
    vec3 b = p1-p;
    vec3 c = p2-p;
    
    float na = length(a);
    float nb = length(b);
    float nc = length(c);
    
    return 2.f * atan(dot(a, cross(b,c)), na*nb*nc + dot(a,b)*nc + dot(b,c)*na + dot(a,c)*nb );

}

float calculateSolidAngle(vec3 x, float levelset, float shift ){
  float angleSum = 0.;
  float dir = 0.f;
  for (int i = 0; i < len2 ;i++){
      
  
    vec3 p0;
    vec3 p1;
    vec3 p2, vX0, vX1, vX2, v, u, c;
    float d, dir;

    /*
    p0 = pts[int(tri[i][0])];
    p1 = pts[int(tri[i][1])];
    p2 = pts[int(tri[i][2])];
    */
    
    p0 = pts[i];
    p1 = pts[(i+1)%(len2)];
    p2 = pts[len2];
    


    float theta =  triangleSA(p0,p1,p2,x);
    // if (theta < 0) theta += 2*M_PI;
    angleSum += theta;
     
  }

  float sa = 2.f*M_PI - angleSum;
  //float solidAng = exteriorAngleSum + (2.f - numPoints) *M_PI;
  return mod(sa - levelset, 4.f*M_PI) + shift ;
}

float map( in vec3 pos )
    
{
    float x = pos[0];
    float y = pos[1];
    float z = pos[2];
    return calculateSolidAngle(pos,0.f, 0.f);
}

//something buggy about this function or like we're not getting accurate harnack results because of it
bool intersectSphere(vec3 ro, vec3 rd, vec3 center, float radius, inout float t0, inout float t1){
  float t = -1.0f;
  vec3 L = ro - center;
  float a = dot(rd,rd);
  float b = 2.0f * dot(rd,L);
  float c = dot(L,L)- dot(radius,radius);
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
        if (t1 < 0.f && t0 < 0.f){
            return false;
        }
       
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

vec3 gradShade( in vec3 p )
{
    vec3 grad = vec3(0.,0.,0.);
    for (int i = 0; i < len-1; i++){
         vec3 p0 = pts[i];
         vec3 p1 = pts[(i+1)%(len-1)];
         vec3 g0 = p0 - p;
         vec3 g1 = p1 - p;
         vec3 n = cross(g1,g0);
         grad += n/dot(n,n) * ((-dot(g0,g1) + dot(g0,g0))/length(g0) + (-dot(g0,g1) + dot(g1,g1))/length(g1));
    }
    return grad/length(grad);
    
}




#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 3
#endif

float getCorrectTime(float t0, float t1){
    if (t0 < 0.f && t1 > 0.f) {
        return t1;
    }
    if (t0  > 0.f) {
        return t0;
    }
}


bool traceVertices(vec3 ro, vec3 rd, inout vec3 finalPos){
    float t0 = 0.0;
    float t1 = 0.0;
    bool didHit;
    for (int i = 0; i < int(numPoints); i ++){
        didHit = intersectSphere(ro, rd, pts[i], .03, t0, t1);
        if (didHit){
            finalPos = ro + getCorrectTime(t0,t1)*rd;
            return true;
        }
     }
    return false;
}

bool raymarch(vec3 ro, vec3 rd, float tmax, inout vec3 finalPos){
    float t0 = 0.0;
    float t1 = 0.0;
    bool didHit;
    for (int i = 0; i < 3; i ++){
        didHit = intersectSphere(ro, rd, pts[i], .05, t0, t1);
        if (didHit){
            finalPos = ro + getCorrectTime(t0,t1)*rd;
            return true;
        }
    }

    
    didHit = rayintersect(ro, rd, t0);
    finalPos = ro + t0*rd;
    return didHit;
        
}



bool closeToVal(float ang, float val){
    float epsilon = .1;
    return ang <= val + epsilon  && ang >= val - epsilon;
}



//For harnack's checks if we're close to our desired level set
bool closeToLevelset(float ang, float levelset){
    vec2 angle = makeComplex(ang);
    vec2 level = makeComplex(levelset);
    float dis = length(angle - level);
    float epsilon = .1;
    return dis <  epsilon;
}



bool raymarchSA(vec3 ro, vec3 rd, float tmax,  inout vec3 finalPos, out bool hitSphere){
    float t = 0.0;
    float t0, t1;
    bool didHit = false;
    
    for( int i=0; i<10000; i++ )
        {
            vec3 pos = ro + t*rd;
            //float h = inPlane(pos);
            
            float h = calculateSolidAngle(pos,0.f,0.f);
            //if( closeToVal(h, 2.*M_PI*cos(iTime)) || t>tmax ) break;
            if( closeToLevelset(h, 2.*M_PI*iTime) || t>tmax   ) break;
            
            t += 0.01;
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
    float shift = 4.f*M_PI;
    float levelset =  2.*M_PI + 1.*M_PI;//*cos(iTime);
    maxSteps = false;
    
  
    tmax = 5.f;
    
    
    pos = ro+t*rd;
    {    
        float loBound = .00001+shift;
        float hiBound = 4.f*M_PI + shift;
        while (t < tmax){
            if (i > 100000){
                maxSteps = true;
                return false;
            }
            float val = calculateSolidAngle(pos,levelset, shift); 
            if (closeToLevelset(val, 0.f + shift)){
                return true;
            }
            
            float R =  closestPointOnLine(pos);
          
            float r = getMaxStep(val, R, loBound, hiBound);
            t += r;
            pos = ro + t*rd;
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
    bool hitSphere = false;
    bool rotate = true;
    
    // camera movement	
	float an;
    if (rotate) an = .8*iTime;
    else an = .8;
	vec3 ro = vec3( 3.0*cos(an), 1.f, 3.0*sin(an) );
	//vec3 ro = vec3( 2., 0.0, 2. );
    vec3 ta = vec3( 0.0, 0.0, 0.0 );
    // camera matrix
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    
    float f1,f2;
    vec3 vertexPos, fnPos;
    
    
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
        
        bool hitVertex  = traceVertices(ro, rd, vertexPos);
       
        if (doHarnacks){
            didHit = harnack(ro,rd, fnPos, maxSteps, iTime);
        } else {
            didHit = raymarchSA(ro,rd,tmax, fnPos, hitSphere );
        }
    
        // shading/lighting	
        vec3 col = vec3(0.0);
        if (maxSteps){
            col =  vec3(1.f,0.0,0.0);
        }
        
        else if (hitVertex && (length(vertexPos-ro) < length(fnPos-ro))){
            vec3 nor = calcNormal(vertexPos);
            float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0));
            col = vec3(0.05,0.5,0.05) + vec3(0.8,0.7,0.5)*dif;
        
        }
        else if( didHit )
        {
            vec3 nor = gradShade(fnPos);
            
            
            float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0));
            if (dot(nor,rd) > 0.){
                col = vec3(0.2,0.3,0.7)*amb + vec3(0.5,0.2,0.9)*dif;
            
            }else{
                col = vec3(0.8,0.3,0.2)*amb + vec3(0.8,0.7,0.2)*dif;
            }
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
