double dt=0.1;
double[] g=new double[]{0,-10}; // gravitational acceleration
int n=100;                     // number of balls
double rho  = 0.1;        // density of balls

double R=250;                   // Half width of the box
double omega=0.;              // angular velocity
double theta=0;                // initial angle
final int viewSize=800;              // size of the viewbox

// Box unitary vectors
double[][] E=new double[4][2];
// Balls parameters
double[] r = new double[n];    // radius
double[] m = new double[n];    // mass
// Balls state and parameters
double[][] x = new double[n][2];  // position
double[][] v = new double[n][2];  // velocity


double norm(double[] v){return Math.sqrt(v[0]*v[0]+v[1]*v[1]);}
double ps(double[] v1,double[] v2){return (v1[0]*v2[0]+v1[1]*v2[1]);}

 
void setup() {
    size(800,800);  // must be equal to viewSize
    // Initialization of the balls
    for(int i=0;i<n;i++){
            r[i]=Math.sqrt(rho)*2*R/sqrt(PI*n);    // radius
      m[i]=PI*r[i]*r[i];            // mass 
            for(int k=0;k<2;k++) x[i][k]=Math.random()*2*(R-r[i])-(R-r[i]);
            v[i][0]=0;v[i][1]=0;
    }
}

void draw() {
  background(1);
  // Let's draw the box and the balls
  translate(viewSize/2,viewSize/2);
  scale(1,-1);
  rotate((float) theta);
  fill(255,255,255);
  rect(-(float) R,-(float) R,2*(float) R,2*(float) R);
  fill(204, 102, 0);
  rotate(-(float) theta); 
  for(int i=0;i<n;i++) circle((float) x[i][0],(float) x[i][1],(float)(2*r[i]));
  
  theta+=omega*dt;
  double[][] E=new double[4][2];
  for(int j=0;j<4;j++) {
    E[j][0]=Math.cos(theta+Math.PI*j/2.);
    E[j][1]=Math.sin(theta+Math.PI*j/2.);
  }
  for(int i=0;i<n;i++){
    // update velocity
    for(int k=0;k<2;k++) v[i][k]+=dt*g[k];
    // update position
    for(int k=0;k<2;k++) x[i][k]+=dt*v[i][k];

    // correction of the velocity and projection wrt the box
    for(int j=0;j<4;j++) 
      if(ps(x[i],E[j])>(R-r[i])){
        // velocity of the wall at x
        double[] v_wall=new double[2];
        v_wall[0]=omega*-x[i][1];v_wall[1]=omega*x[i][0];
        // relative velocity
        double[] v_relative=new double[2];
        for(int k=0;k<2;k++) v_relative[k]=v[i][k]-v_wall[k];
        // projection de la vitesse relative
        double v_normal=ps(v_relative,E[j]);
        if(v_normal>0) 
          for(int k=0;k<2;k++) {
            v_relative[k]-=2*v_normal*E[j][k];
            v[i][k]=v_wall[k]+v_relative[k];
          }
        for(int k=0;k<2;k++) x[i][k]-=(ps(x[i],E[j])-(R-r[i]))*E[j][k];
      }
  }
  // collision between balls
  // Note : Probably better to make a copy of v[] before... 
  //    ... it would be cleaner (but it is not clear it would lead to noticeable improvements).
  for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
  double[] dx={x[i][0]-x[j][0],x[i][1]-x[j][1]};// relative position of the balls
  double[] dv={v[i][0]-v[j][0],v[i][1]-v[j][1]};// relative velocity of the ball
  double norm_dx=norm(dx);
  if(norm_dx<r[i]+r[j])  // balls overlapping
  if(ps(dx,dv)<0)      // ball colliding
    {
      /*
      m1 V1=m1v1 + P
      m2 V2=m2v2 - P
      m1 V1^2 + m2 V2^2 = m1v1^2 + m2v2^2 +2P(v1-v2) + (1/m1+1/m2)P^2
      2(v1-v2) = I-(m1+m2)/(m1m2) P
      P=2m1m2/(m1+m2)*(v1-v2)
      */
      double[] n=new double[]{dx[0]/norm_dx,dx[1]/norm_dx};
      double P=-2*m[i]*m[j]/(m[i]+m[j])*ps(dv,n);
      for(int k=0;k<2;k++){
        v[i][k]+=P*n[k]/m[i];
        v[j][k]-=P*n[k]/m[j];
      }
    }
  }
}
