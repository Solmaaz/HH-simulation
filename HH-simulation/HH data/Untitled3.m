 %i n i t i a l i z e pa r ame t e r s
 dt=0.5;
 d =8;
 a =0.02;
 c =-65; 
 b =0.2;

 %r e s e r v e memory
 T =ceil(1000/dt) ;
 v =zeros(T,1) ;
 u =zeros(T,1) ;
 v(1)=-70;
 u(1)= -14;

 nin=100 ; %NEW
 rate=2*1e-3;%to [ms ] %NEW
 taug=10 ; %NEW
gin=zeros(nin,1) ; %NEW 
Ein=zeros(nin,1) ; %NEW
win=0.07*ones(1,nin) ; %NEW

 %f o r?l o o p o v e r t ime
 for t=1:T-1;
 %Po i s s o n i n p u t
 if t*dt>200 && t*dt<700
 p= rand(nin,1)<rate*dt ; %NEW
 else
 p=0 ; %NEW
 end
 %conduc tanc e update o f g i n
 gin=gin+p ;
 Iapp=win*(gin.*Ein) ;
 Iapp=Iapp-(win*gin).*v(t) ;
 gin=(1-dt/taug)*gin ;

 if v(t)<35
 %update ODE
 dv=(0.04*v(t)+5)*v(t)+140-u(t);
 v(t+1)=v(t)+(dv+Iapp)*dt;
 du=a*(b*v(t)-u(t));
 u(t+1)=u(t)+dt*du;
 else
 %s p i k e s
 v(t)=35;
 v(t+1)=c;
 u(t+1)=u(t)+d;
 end

 end

 %p l o t v o l t a g e t r a c e
 plot((0:T-1)*dt,v,'b') ;
 xlabel('Time[ms]') ;
 ylabel( 'Membrane voltage [mV]' ) ;