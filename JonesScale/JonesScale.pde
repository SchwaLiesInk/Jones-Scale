/**
 * The Jone Set
 * by Gerwyn Jones.  
 * 
 * Simple rendering of the Mandelbrot set.
 */

import java.awt.event.*;
// Establish a range of values on the complex plane
// A different range will allow us to "zoom" in or out on the fractal
// float xmin = -1.5; float ymin = -.1; float wh = 0.15;

double frame=1;
double speed=0.0001;
double step=speed;
double anim=0.0001;
double zoom=5.0;
double zx=0;
double zy=0;
double xz=0;
double yz=0;
double xm=0;
double ym=0;
double setc=1.0/3.0;
double set=0;
boolean released=false;
PImage mset;

// Make sure we can write to the pixels[] array.
// Only need to do this once since we don't do any other drawing.

void mouseReleased() {
  double dx=(ym/height);
  double dy=(xm/width); 
  yz=dy-0.5;xz=dx-0.5;
  
}
void mouseMoved() {
  xm=mouseX;
  ym=mouseY;
}
double t=1;
double tz=1;
double []formula(double za,double zb,double ca,double cb) {
  //
  double radius=Math.sqrt(ca*ca+cb*cb);
  double ir=1.0/radius;
  double vol=4.0/3.0*radius*radius*radius;
  double rx=1.0/(radius/t);
  double ry=(1.0/t-ir);
  double na=za*rx;
  double nb=zb*rx;
  na=na/(((na/za-za*ry)*ir)*vol);
  nb=nb/(((nb/zb-zb*ry)*ir)*vol);
  //
  na=(na+ca+tz);
  nb=(nb+cb+tz);
  return new double[]{na,nb};
}
void draw() {
  if (keyCode==UP) {
    zoom*=0.9;
  } else if (keyCode==DOWN) {
    zoom*=1.1;
  }
  
  if (keyCode==LEFT) {
    anim=-step*Math.sqrt(zoom);
  } else if (keyCode==RIGHT) {
    anim=step*Math.sqrt(zoom);
  }
  if(keyCode==java.awt.event.KeyEvent.VK_W){
    t*=1.1;
  }else
  if(keyCode==java.awt.event.KeyEvent.VK_S){
    t*=0.9;
  }
  
  if(keyCode==java.awt.event.KeyEvent.VK_E){
    tz+=0.1;
  }else
  if(keyCode==java.awt.event.KeyEvent.VK_D){
    tz-=0.1;
  }
  frame+=anim;
  
  if(frame>1){step=-speed;frame=1;}else
  if(frame<-1){step=speed;frame=-1;}
   
  MSet();
  yz=0;xz=0;
  mset.updatePixels();
  image(mset, 0, 0);
}
void setup() {

  size(800, 800);
  background(0);
  mset=new PImage(width, height);
  MSet();
  mset.updatePixels();
}
void MSet() {;
  zx+=xz*zoom;
  zy+=yz*zoom;
  
    //t*=1.1;
    //zoom*=1.1;
  double xmin = ((-zoom+zx)+zoom*0.5);
  double ymin = ((-zoom+zy)+zoom*0.5);
  double w = (zoom);
  double h = (zoom);
  // Maximum number of iterations for each point on the complex plane
  int maxiterations = 16;

  // x goes from xmin to xmax
  double xmax = xmin + w;
  // y goes from ymin to ymax
  double ymax = ymin + h;
  // Calculate amount we increment x,y for each pixel
  double dx =(xmax - xmin) / (width);
  double dy =(ymax - ymin) / (height);

  // Start y
  double y = ymin;
  for (int j = 0; j < width; j++) {
    //int jc=j-(width>>1);
    // Start x
    double x = xmin;
    for (int i = 0; i < height; i++) {
    //int ic=i-(height>>1);
      //double r=0.456789/(Math.sqrt(x*x+y*y)*0.5);
      //double r=1.0/(Math.sqrt(x*x+y*y)*0.5);
      //int alpha=128-(int)(Math.sqrt(ic*ic+jc*jc)*128)/width;
      //if(alpha<0){alpha=0;}
      // Now we test, as we iterate z = z^2 + cm does z tend towards infinity?
      double ax=x;
      double ay=y;
      //double tr=1.0;//Math.log10(t)/(x*x+y*y);//(x*x+y*y)/t;
      double bx=(ax*Math.cos(0.7853981634)-ay*Math.sin(0.7853981634));
      double by=(ay*Math.cos(0.7853981634)+ax*Math.sin(0.7853981634));
      double a = bx;//Math.abs(x);//(x*frame);//*r;
      double b = by;//Math.abs(y);//(y*frame);//*r;
      int n = 0;
      while (n < maxiterations) {
        
        //
        double na[]=formula(a,b,bx,by);
        a =na[0];
        b =na[1];
        
        // Infinty in our finite world is simple, let's just consider it 16
        if (na[0] +na[1] > maxiterations) {
         break;  // Bail
        }
        n++;
      }

      // We color each pixel based on how long it takes to get to infinity
      // If we never got there, let's pick the color black
      //if (n == maxiterations) {
        //mset.pixels[j+i*width] = color(0,alpha);
      //} else {
        // Gosh, we could make fancy colors here if we wanted
        /*double t10=Math.abs(Math.log10(t)/Math.log(t));
        double r=Math.sqrt(a*a + b*b);
        r=r/Math.pow(t,2.718);
        double r10=Math.log10(r);
        r10=Math.abs(Math.pow(r10,2))/t10;
        int rd=(int)((1+Math.sin(r))*r10)%255;
        int gr=(int)((1+Math.cos(r))*r10)%255;
        int bl=(int)((1-Math.sin(r))*r10)%255;*/
        double r=Math.sqrt(a*a + b*b)/t;
        int rd=(int)((1+Math.sin(r))*128);
        int gr=(int)((1+Math.cos(r))*128);
        int bl=(int)((1-Math.sin(r))*128);
        mset.pixels[j+i*width] = color( rd,gr,bl,255);
       //mset.pixels[j+i*width] = color( rd,gr,bl,64);
       //mset.pixels[j+i*width] = color( rd,gr,bl,64);
      //}
      x += dx;
    }
    y += dy;
  }
  set+=setc;
  if(set>255){
  set-=255;
  }//else if(set==0){
  //setc=-setc;
  //}
}