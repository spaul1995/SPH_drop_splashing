float g = 9.8;
float dt = 0.05;
float k = 100;
float limit = 0.1;
float sigma = -1;
float t = 0;
float viscosity = 10;
float densitykernel = 7;
float pressurekernel = 1.5*densitykernel;
float viscositykernel = 1.5*densitykernel;
float colorkernel = 1.001*densitykernel;
float tensionkernel = 2*densitykernel;
float initialspacing = 1*densitykernel;
float boxclearance = 50;
float boxclearancey = 100;
float radius = 6;
float radius_b = 2;
float screen_l = 500;
int nos =625;
int nos_l = (int)Math.sqrt(nos);
PVector firstball = new PVector(screen_l/2-0.5*(nos_l-1)*initialspacing,screen_l/2-0.5*(nos_l-1)*1*initialspacing);
float boundaryballseperation = 0.5*initialspacing;
//int nos_b_l = (int)((screen_l-firstball.y-boxclearance)/boundaryballseperation)+1;
int nos_b_b = (int)((screen_l/2)/boundaryballseperation);
int nos_b_l = (int)(pow((pow((screen_l-firstball.y-boxclearance),2)+pow((screen_l-nos_b_b*boundaryballseperation-2*boxclearancey)/2,2)),0.5)/boundaryballseperation)+1;
int nos_b_r = nos_b_l-1;
//int nos_b_b = (int)((screen_l-2*boxclearancey)/boundaryballseperation);

//int nosbtotal=nos_b_l+nos_b_b+nos_b_r;
int nosbtotal=nos_b_l+nos_b_r+nos_b_b;
float[][] colorf = new float[(int)screen_l][(int)screen_l];
float[][] gradientcolorfx = new float[(int)screen_l][(int)screen_l];
float[][] gradientcolorfy = new float[(int)screen_l][(int)screen_l];
float[][] laplaciancolorf = new float[(int)screen_l][(int)screen_l];
PVector[] r = new PVector[nos+nosbtotal];
PVector[] bracbillf = new PVector[nos];
float[] newsurfacetensionforcex = new float[nos+nosbtotal];
float[] newsurfacetensionforcey = new float[nos+nosbtotal];
float[] rho = new float[nos+nosbtotal];
float[] pressure = new float[nos+nosbtotal];
float[] grad_pressurex = new float[nos+nosbtotal];
float[] grad_pressurey = new float[nos+nosbtotal];
float[] laplacian_ux = new float[nos+nosbtotal];
float[] laplacian_uy = new float[nos+nosbtotal];
float[] m = new float[nos+nosbtotal];
PVector[] v = new PVector[nos+nosbtotal];
float gridsize = densitykernel;
float soundsq = 0.005;
float rho_initial = 2;
float m_initial = rho_initial*3.14*densitykernel*densitykernel;
int gridnos = (int)(screen_l/gridsize);
int maxparticlespergrid = 500;
int[][][] gridstore = new int[gridnos][gridnos][maxparticlespergrid];
int[][] newcounter= new int [gridnos][gridnos];
Table table;
float sin= (screen_l-firstball.y-boxclearance)/pow((pow((screen_l-firstball.y-boxclearance),2)+pow(((screen_l-nos_b_b*boundaryballseperation-2*boxclearancey)/2),2)),0.5);
float cos= ((screen_l-nos_b_b*boundaryballseperation-2*boxclearancey)/2)/pow((pow((screen_l-firstball.y-boxclearance),2)+pow(((screen_l-nos_b_b*boundaryballseperation-2*boxclearancey)/2),2)),0.5); 




void setup() {
  size(500,500);
  table = new Table();
  table.addColumn("time");
  table.addColumn("bracbill");
  //table.addColumn("message");
  for (int j=0;j<500;j++){
      for (int i=0;i<500;i++){
        colorf[j][i]=0;
      }
    }

  // creating the initial particles
    for (int j=0;j<nos_l;j++){
      for (int i=0;i<nos_l;i++){
      r[(nos_l*j+i)]=new PVector(firstball.x+(i)*initialspacing,firstball.y+(j)*initialspacing);
      v[nos_l*j+i]=new PVector(0,0);
      bracbillf[nos_l*j+i]=new PVector(0,0);
      rho[(nos_l*j+i)] = rho_initial;
      m[(nos_l*j+i)] = m_initial;
    }
    }
  // creating the boundary walls
    PVector firstballofleftboundary = new PVector(boxclearancey,firstball.y);
    //PVector firstballofrightboundary = new PVector(screen_l-boxclearance,firstball.y);
    int counter=nos;
    int counter1 = 0;
    //while (counter1<nos_b_l){
    //  r[counter]=new PVector(firstballofleftboundary.x,firstballofleftboundary.y+counter1*boundaryballseperation);
    //  v[counter]=new PVector(0,0);
    //  rho[counter] = rho_initial;
    //  m[counter] = m_initial;
    //  counter+=1;
    //  counter1+=1;
    //}
    //while (counter1<nos_b_l+nos_b_b){
    //  r[counter]=new PVector(firstballofleftboundary.x+(counter1-nos_b_l+1)*boundaryballseperation,screen_l-boxclearance);
    //  v[counter]=new PVector(0,0);
    //  rho[counter] = rho_initial;
    //  m[counter] = m_initial;
    //  counter+=1;
    //  counter1+=1;
    //}   
    //while (counter1<nos_b_l+nos_b_b+nos_b_r){
    //  r[counter]=new PVector(firstballofleftboundary.x+(nos_b_b)*boundaryballseperation,screen_l-boxclearance-(counter1-nos_b_l-nos_b_b+1)*boundaryballseperation);
    //  v[counter]=new PVector(0,0);
    //  rho[counter] = rho_initial;
    //  m[counter] = m_initial;
    //  counter+=1;
    //  counter1+=1;
   
    while (counter1<nos_b_l){
      r[counter]=new PVector(firstballofleftboundary.x+counter1*boundaryballseperation*cos,firstballofleftboundary.y+counter1*boundaryballseperation*sin);
      //print(r[counter],"\n");
      v[counter]=new PVector(0,0);
      rho[counter] = rho_initial;
      m[counter] = m_initial;
      counter+=1;
      counter1+=1;
    }
    while (counter1<nos_b_l+nos_b_b){
      r[counter]=PVector.add(r[nos+nos_b_l-1],new PVector((counter1-nos_b_l)*boundaryballseperation,0));
      v[counter]=new PVector(0,0);
      rho[counter] = rho_initial;
      m[counter] = m_initial;
      counter+=1;
      counter1+=1;
    }    
    //counter1+=1;
    while (counter1<nosbtotal){
      r[counter]=PVector.add(r[nos+nos_b_l+nos_b_b-1],new PVector((counter1-nos_b_l-nos_b_b)*boundaryballseperation*cos,-(counter1-nos_b_l-nos_b_b)*boundaryballseperation*sin));
      print(r[nos+nos_b_l+nos_b_b-1],"\n");
      //print(r[counter],"\n");
      v[counter]=new PVector(0,0);
      rho[counter] = rho_initial;
      m[counter] = m_initial;
      counter+=1;
      counter1+=1;
    }    
fill(color(30,144,255));
  stroke(0);
  for (int j=0;j<nos+nosbtotal;j++){
  ellipse(r[j].x,r[j].y,radius,radius);
  }      
      
      
      
}



void draw() {
  background(255); 
  gridassigning();
  es_density();
  es_tension();
  es_surfacetension();
  es_pressure();
  es_gradpressure();
  es_viscosity();
  move();
  display();
  saveFrame("images/squares-######.png");
}


void gridassigning() {
  gridstore=new int[gridnos][gridnos][maxparticlespergrid];
  newcounter= new int [gridnos][gridnos];
  for (int j=0;j<nos+nosbtotal;j++){
    gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)]]=j;
    newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)]++;
    if (newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)]>maxparticlespergrid){
      print("out");
    }
  }
  //print(newcounter[4][4], "    ");
}

void adddensity(int a, int b){
  rho[a]+=m[b]*3/(3.14*pow(densitykernel,6))*(pow((densitykernel*densitykernel-PVector.sub(r[a],r[b]).magSq()),2));
  }
void addtension(int a, int b){
  float kernelte=0;
  if (PVector.sub(r[a],r[b]).mag()<tensionkernel/3){
    kernelte=(63/(478*3.14*tensionkernel*tensionkernel))*(66-pow((3-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5)-6*pow((2-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5)+15*pow((2-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5));
  }
  else if (PVector.sub(r[a],r[b]).mag()>tensionkernel/3 && PVector.sub(r[a],r[b]).mag()<2*tensionkernel/3){
    kernelte=(63/(478*3.14*tensionkernel*tensionkernel))*(66-pow((3-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5)-6*pow((2-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5));
  }
  else if (PVector.sub(r[a],r[b]).mag()>2*tensionkernel/3){
    kernelte=(63/(478*3.14*tensionkernel*tensionkernel))*(66-pow((3-3*(PVector.sub(r[a],r[b]).mag()/tensionkernel)),5));
  }
  newsurfacetensionforcex[a]+=PVector.mult(PVector.sub(r[a],r[b]),(sigma*kernelte)).x;
  newsurfacetensionforcey[a]+=PVector.mult(PVector.sub(r[a],r[b]),(sigma*kernelte)).y;  
}
void addgradpx(int a, int b){
  grad_pressurex[a]+=(m[b]/rho_initial)*(pressure[a]+pressure[b])/2*(-30/(3.14*pow(pressurekernel,5)))*pow((pressurekernel-PVector.sub(r[a],r[b]).mag()),2)*(r[b].x-r[a].x)/(PVector.sub(r[a],r[b]).mag());
}
void addgradpy(int a, int b){
    grad_pressurey[a]+=(m[b]/rho_initial)*(pressure[a]+pressure[b])/2*(-30/(3.14*pow(pressurekernel,5)))*pow((pressurekernel-PVector.sub(r[a],r[b]).mag()),2)*(r[b].y-r[a].y)/(PVector.sub(r[a],r[b]).mag());
}
void addlaplacianux(int a, int b){
 laplacian_ux[a]+=viscosity*(m[b]/rho_initial)*(v[b].x-v[a].x)*(45/(3.14*pow(viscositykernel,5)))*(viscositykernel-PVector.sub(r[a],r[b]).mag());
}
void addlaplacianuy(int a, int b){
  //print(v[a].y);
  laplacian_uy[a]+=viscosity*(m[b]/rho_initial)*(v[b].y-v[a].y)*(45/(3.14*pow(viscositykernel,5)))*(viscositykernel-PVector.sub(r[a],r[b]).mag());
}
void addcolor(int a, int c, int b){
  colorf[a][c]+=(m[b]/rho_initial)*3/(3.14*pow(colorkernel,6))*(pow((colorkernel*colorkernel-PVector.sub(new PVector(a,c),r[b]).magSq()),2));
  laplaciancolorf[a][c]+=(m[b]/rho_initial)*(-6/(3.14*pow(colorkernel,6)))*(2*colorkernel*colorkernel-4*PVector.sub(new PVector(a,c),r[b]).magSq());
  gradientcolorfx[a][c]+=(m[b]/rho_initial)*(-6/(3.14*pow(colorkernel,6)))*PVector.sub(new PVector(a,c),r[b]).mag()*(colorkernel*colorkernel-PVector.sub(new PVector(a,c),r[b]).magSq())*(PVector.sub(new PVector(a,c),r[b]).normalize().x);
  gradientcolorfy[a][c]+=(m[b]/rho_initial)*(-6/(3.14*pow(colorkernel,6)))*PVector.sub(new PVector(a,c),r[b]).mag()*(colorkernel*colorkernel-PVector.sub(new PVector(a,c),r[b]).magSq())*(PVector.sub(new PVector(a,c),r[b]).normalize().y);
  ////print(rho[a], "   ");
}

void es_pressure() {
  for (int j=0;j<nos+nosbtotal;j++){
      //Tait equation
      pressure[j]=rho_initial*soundsq*(pow((rho[j]/rho_initial),7)-1);
      //print( (rho[j]/rho_initial), "   " , pressure[j], "   ");
  }
}



void colorfield() {
for (int m=50;m<450;m++){
      for (int n=50;n<450;n++){
        
        for (int i=0;i<newcounter[(int)(m/gridsize)+1][(int)(n/gridsize)];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)][i]);
         //print("east    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)-1][(int)(n/gridsize)];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)][i]);
         //print("west    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)][(int)(n/gridsize)-1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)][(int)(n/gridsize)-1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)][(int)(n/gridsize)-1][i]);
         //print("north    ");
      }
    }    
        for (int i=0;i<newcounter[(int)(m/gridsize)][(int)(n/gridsize)+1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)][(int)(n/gridsize)+1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)][(int)(n/gridsize)+1][i]);
         //print("south    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)+1][(int)(n/gridsize)-1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)-1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)-1][i]);
         //print("northeast    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)+1][(int)(n/gridsize)+1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)+1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)+1][(int)(n/gridsize)+1][i]);
         //print("southeast    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)-1][(int)(n/gridsize)-1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)-1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)-1][i]);
         //print("northwest    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)-1][(int)(n/gridsize)+1];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)+1][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)-1][(int)(n/gridsize)+1][i]);
         //print("southwest    ");
      }
    }
        for (int i=0;i<newcounter[(int)(m/gridsize)][(int)(n/gridsize)];i++){
      if (PVector.sub(new PVector(m,n),r[gridstore[(int)(m/gridsize)][(int)(n/gridsize)][i]]).mag()<colorkernel){
         addcolor(m, n, gridstore[(int)(m/gridsize)][(int)(n/gridsize)][i]);
         //print("center    ");
      }
    }    
      }
}

//calc_gradcolour();
//calc_lapcolour();
//showcolorfield();

//print(colorf[(int)r[0].x][(int)r[0].y],"\t",colorf[(int)r[nos-1].x][(int)r[nos-1].y], "\n");
//print(new PVector(gradientcolorfx[(int)r[0].x][(int)r[0].y],gradientcolorfy[(int)r[0].x][(int)r[0].y] ),"\t",new PVector((gradientcolorfx[(int)r[nos-1].x][(int)r[nos-1].y]), (gradientcolorfy[(int)r[nos-1].x][(int)r[nos-1].y] )), "\n");
//print(laplaciancolorf[(int)r[0].x][(int)r[0].y],"\t",laplaciancolorf[(int)r[nos-1].x][(int)r[nos-1].y], "\n");
}
//void calc_gradcolour() {
//  for (int m=50;m<450;m++){
//      for (int n=50;n<450;n++){
//        //gradientcolorfx[m][n]=(-colorf[m+2][n]+8*colorf[m+1][n]-8*colorf[m+1][n]+colorf[m-2][n])/12;
//        //gradientcolorfy[m][n]=(-colorf[m][n+2]+8*colorf[m][n+1]-8*colorf[m][n+1]+colorf[m][n-2])/12;
//        gradientcolorfx[m][n]=(colorf[m+2][n]-2*colorf[m][n]+colorf[m-2][n])/4;
//        gradientcolorfy[m][n]=(colorf[m][n+2]-2*colorf[m][n]+colorf[m][n-2])/4;
//      }
//  }
//}

//void calc_lapcolour() {
//  for (int m=50;m<450;m++){
//      for (int n=50;n<450;n++){
//        laplaciancolorf[m][n]=(gradientcolorfx[m+1][n]-gradientcolorfx[m-1][n])/2+(gradientcolorfx[m][n-1]-gradientcolorfx[m][n+1])/2;
//                if (t==0 && m==250){
//    //print(new PVector(gradientcolorfx[250][n],gradientcolorfy[250][n]).mag(),"\n");
//    //print(colorf[250][n],"\n");
//    }
//      }
//  }
//}
void showcolorfield() {
  for (int m=50;m<450;m++){
      for (int n=50;n<450;n++){
        if (colorf[m][n]!=0){
set(m, n, int(255));

//print(colorf[m][n],"    ");
        }

      }
      }
}


void es_surfacetension() {
  
  for (int j=0;j<nos;j++){
    //if (t==0){
    //print(new PVector(gradientcolorfx[(int)(r[j].x)][(int)(r[j].y)],gradientcolorfy[(int)(r[j].x)][(int)(r[j].y)]).mag(),"\n");
    //}
    //print(colorf[(int)(r[j].x)][(int)(r[j].y)],"\n");
    if (new PVector(gradientcolorfx[(int)(r[j].x)][(int)(r[j].y)],gradientcolorfy[(int)(r[j].x)][(int)(r[j].y)]).mag()>limit){
    bracbillf[j] = PVector.mult(new PVector(gradientcolorfx[(int)(r[j].x)][(int)(r[j].y)],gradientcolorfy[(int)(r[j].x)][(int)(r[j].y)]).normalize(),
    -sigma*laplaciancolorf[(int)(r[j].x)][(int)(r[j].y)]);
    //print(j,"\t",new PVector(gradientcolorfx[(int)(r[j].x)][(int)(r[j].y)],gradientcolorfy[(int)(r[j].x)][(int)(r[j].y)]),"\n");
    //print(j,"\t",bracbillf[j],"\n");
    }

  }
      //print(1,"\t",bracbillf[0],"\t",nos-1,"\t",bracbillf[nos-1],"\n");
  //TableRow newRow = table.addRow();
  //newRow.setFloat("time", t);
  //newRow.setFloat("bracbill", bracbillf[1].x);
  //saveTable(table, "data/new.csv");
  t+=dt;
}





void es_viscosity() {
for (int j=0;j<nos+nosbtotal;j++){
    //adding density for east box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
      }
    }
        //adding density for west box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
      }
    }
        //adding density for north box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
      }
    }
        //adding density for south box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
      }
    }
        //adding density for northeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
      }
    }
        //adding density for southeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
      }
    }
        //adding density for southwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
      }
    }
        //adding density for northwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
      }
    }
        //adding density for center box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]]).mag()<viscositykernel){
         addlaplacianux(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
         addlaplacianuy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
      }
    }
    //print(laplacian_uy[nos_l*(nos_l-1)+nos_l/2],"    ", laplacian_uy[nos_l*(nos_l-1)-nos_l/2],"\n");
  }  
  
}





void move() {
  for (int j=0;j<nos;j++){
    //v[j].x+=dt*(viscosity*laplacian_ux[j]-1/(rho[j])*grad_pressurex[j]+bracbillf[j].x);
    //v[j].y+=dt*(viscosity*laplacian_uy[j]-1/(rho[j])*grad_pressurey[j]+g+bracbillf[j].y);
    
    //v[j].x+=dt*(bracbillf[j].x);
    //v[j].y+=dt*(bracbillf[j].y);    
    
    v[j].x+=dt*(+1/(rho[j])*grad_pressurex[j]+laplacian_ux[j]+newsurfacetensionforcex[j]);
    v[j].y+=dt*(+1/(rho[j])*grad_pressurey[j]+g+laplacian_uy[j]+newsurfacetensionforcey[j]);
    
    
    //bracbillf[j].x=0;
    //bracbillf[j].y=0;
    //r[j].x+=dt*v[j].x;
    //r[j].y+=dt*v[j].y;
    //r[j].x+=dt*1;
    r[j].x+=dt*v[j].x;
    r[j].y+=dt*v[j].y;
    //print(rho[j], "   ");
  }
    print(newsurfacetensionforcex[0], "\n");
    print("velocity", "\t",v[0].x,"\n");
  //print(v[nos_l*(nos_l-1)+nos_l/2].y,"\n");
rho = new float[nos+nosbtotal];
for (int j=0;j<nos+nosbtotal;j++){
  rho[j]=rho_initial;
}
pressure = new float[nos+nosbtotal];
grad_pressurex = new float[nos+nosbtotal];
grad_pressurey = new float[nos+nosbtotal];
laplacian_ux = new float[nos+nosbtotal];
laplacian_uy = new float[nos+nosbtotal];
colorf = new float[500][500];
laplaciancolorf = new float[500][500];
gradientcolorfx = new float[500][500];
gradientcolorfy = new float[500][500];
newsurfacetensionforcex = new float[nos+nosbtotal];
newsurfacetensionforcey = new float[nos+nosbtotal];
}





void display() {
  fill(color(30,144,255));
  stroke(0);
  for (int j=0;j<nos+nosbtotal;j++){
  ellipse(r[j].x,r[j].y,radius,radius);
  }
//  loadPixels();
//for (int n=50;n<450;n++){
//      for (int m=50;m<450;m++){
//        if (colorf[m][n]!=0){
//pixels[400*n+m] = (int)(255*(1-colorf[m][n]));
////print(colorf[m][n],"    ");
//        }
//      }
//      }
//updatePixels();
   
  //for (int j=0;j<gridnos;j++){
  //   line(j*gridsize,0,j*gridsize,screen_l);
  //   line(0,j*gridsize,screen_l,j*gridsize);
  //}
}


void es_gradpressure() {
  for (int j=0;j<nos+nosbtotal;j++){
    grad_pressurex[j]=0;
    //adding density for east box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
         //print("oops1     ");
      }
    }
        //adding density for west box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
         //print("oops1     ");
      }
    }
        //adding density for north box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
         //print("oops1     ");
      }
    }
        //adding density for south box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
         //print("oops1     ");
      }
    }
        //adding density for northeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops1     ");
      }
    }
        //adding density for southeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops1     ");
      }
    }
        //adding density for southwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops1     ");
      }
    }
        //adding density for northwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]]).mag()<pressurekernel){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops1     ");
      }
    }
        //adding density for center box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]]).mag()<pressurekernel){
        if (j!=gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]){
         addgradpx(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
         addgradpy(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
         //print("oops1     ");
        }
      }
    }
    //print(grad_pressurex[j],"    ");
    //print(grad_pressurey[j],"    ");
  }
  //print(pressure[884],"    ", grad_pressurey[884],"    ", pressure[854],"    ", grad_pressurey[854],"\n");
}

void es_density() {
  for (int j=0;j<nos+nosbtotal;j++){
    //adding density for east box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");
         //print("east    ");
      }
    }
        //adding density for west box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");         
         //print("west    ");
      }
    }
        //adding density for north box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for south box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for northeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for southeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for southwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for northwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]]).mag()<densitykernel){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for center box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]]).mag()<densitykernel){
        if (j!=gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]){
         adddensity(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");
       }
      }
    }
    //print(rho[j], "   ");  
  }
}
void es_tension() {
  for (int j=0;j<nos;j++){
    //adding density for east box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");
         //print("east    ");
      }
    }
        //adding density for west box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");         
         //print("west    ");
      }
    }
        //adding density for north box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for south box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for northeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for southeast box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)+1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for southwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)+1][i]);
         //print("oops     ");
      }
    }
        //adding density for northwest box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]]).mag()<tensionkernel){
         addtension(j, gridstore[(int)(r[j].x/gridsize)-1][(int)(r[j].y/gridsize)-1][i]);
         //print("oops     ");
      }
    }
        //adding density for center box
    for (int i=0;i<newcounter[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)];i++){
      if (PVector.sub(r[j],r[gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]]).mag()<tensionkernel){
        if (j!=gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]){
         addtension(j, gridstore[(int)(r[j].x/gridsize)][(int)(r[j].y/gridsize)][i]);
         //print("oops     ");
       }
      }
    }
    //print(rho[j], "   ");  
  }
}
