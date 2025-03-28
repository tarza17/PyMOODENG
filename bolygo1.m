%% Bolygomozgas
% A Fold Nap koruli keringeset szimulalom kulonbozo modszerekkel. A Fold palyaja nem annyira erdekes, mert kozel kor alaku. A kezdeti feltetelek allitasaval lathatjuk, hogy a modszerek kiszamoljak az ellipszis palyakat is.
%
% A dinamika alaptorvenye:
%
% $\mathbf{F} = \dot{\mathbf{p}} = m\mathbf{a} = m\ddot{\mathbf{r}}$
%
% ahol egy $M$ tömegû origóban rögzített test által az $m$ tömegû testre ható erõ
%
% $\mathbf{F} = -G\frac{mM}{r^2}\frac{\mathbf{r}}{r}$
%
% ahol $G$ a gravitációs állandó. Tehát a test mozgásegyenlete:
%
% $\ddot{\mathbf{r}} = -G\frac{m}{r^3}\mathbf{r}$
%
% Ekkor bevezetve a $\mathbf{v}=\dot{\mathbf{r}}$ sebesseget az alabbi ketismeretlenes differencialegyenlet-rendszert kapjuk
%
% $\dot{\mathbf{v}} = -G\frac{M}{r^3}\mathbf{r}$
%
% $\dot{\mathbf{r}} = \mathbf{v}$
%
% Az egyes modszereknel beallithato, hogy hany evig szimulalja a mozgast a
% "year" valtozoval illetve a lepeskozt a "day" valtozoval (utobbi nem allithato az
% ode45-nel).


%% Euler modszer
%

clc;
clear all;

%Allandok
m=5.972e+24;        %Fold tomege            (kg)
M=1.989e+30;        %Nap tomege             (kg)
G=6.67408e-11;      %Gravitacios allando    (m^3*kg^-1*s^-2)

%Kezdeti feltetelek
vx(1)=0;            %Napkozelben a sebesseg (m/s)
vy(1)=30.29e+3;     %Napkozelben a sebesseg (m/s)

rx(1)=147.1e+9;     %Napkozelben a tavolsag (m)
ry(1)=0;            %Napkozelben a tavolsag (m)

year=5;             %Hany evig keringjen a Fold a Nap korul
day=1;              %Lepeskoz, hany naponta szamoljon

t=60*60*24*365*year;%Keringesi ido          (secundum)
h=60*60*24*day;     %Lepeskoz               (secundum)
time(1)=0;          %Kezdettol eltelt ido   (secundum)

%Szamitas Euler modszer

ax(1)=-G*M*rx(1)/(rx(1)^2+ry(1)^2)^1.5;     %Sebesseg   (m*s^-2)
ay(1)=-G*M*ry(1)/(rx(1)^2+ry(1)^2)^1.5;     %Sebesseg   (m*s^-2)

E_mozg(1)=(vx(1)^2+vy(1)^2)/2*m;            %Mozgasi energia
E_hely(1)=-G*M*m/(rx(1)^2+ry(1)^2)^0.5;     %Helyzeti energia
E_ossz(1)=E_mozg(1)+E_hely(1);              %Osszenergia

L(1)=m*norm(cross([rx(1),ry(1),0],[vx(1),vy(1),0]));    %Perdulet

for i = 2:(t/h)
    ax(i)=-G*M*rx(i-1)/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    ay(i)=-G*M*ry(i-1)/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    
    vx(i)=vx(i-1)+ax(i)*h;
    vy(i)=vy(i-1)+ay(i)*h;
    
    rx(i)=rx(i-1)+vx(i)*h;
    ry(i)=ry(i-1)+vy(i)*h;
    
    E_mozg(i)=(vx(i)^2+vy(i)^2)/2*m;
    E_hely(i)=-G*M*m/(rx(i)^2+ry(i)^2)^0.5;
    E_ossz(i)=E_mozg(i)+E_hely(i);
    
    L(i)=m*norm(cross([rx(i), ry(i), 0],[vx(i), vy(i), 0]));
    
    time(i)=time(i-1)+h;
end

%Kiiratas

figure;
str=sprintf('Euler modszer \n(keringesi ido %d ev, lepeskoz %d nap)', year, day);
sgtitle(str,'Interpreter', 'latex');
subplot(331);
hold on;
%axis equal;
plot(rx, ry, 'b');
plot(0, 0, 'oy', 'MarkerFaceColor', 'y');
xlabel('$r_x(t)$', 'Interpreter', 'latex');
ylabel('$r_y(t)$', 'Interpreter', 'latex');
title('Palyagorbe', 'Interpreter', 'latex');

subplot(332);
hold on;
plot(time, vx, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_x(t)$', 'Interpreter', 'latex');
title('X sebesseg-ido fv', 'Interpreter', 'latex');

subplot(333);
hold on;
plot(time, vy, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_y(t)$', 'Interpreter', 'latex');
title('Y sebesseg-ido fv', 'Interpreter', 'latex');

subplot(334);
hold on;
plot(time, E_mozg, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{mozg}(t)$', 'Interpreter', 'latex');
title('Mozgasi energia', 'Interpreter', 'latex');

subplot(335);
hold on;
plot(time, E_hely, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{hely}(t)$', 'Interpreter', 'latex');
title('Helyzeti energia', 'Interpreter', 'latex');

subplot(336);
hold on;
plot(time, E_ossz, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{ossz}(t)$', 'Interpreter', 'latex');
title('Osszenergia', 'Interpreter', 'latex');

subplot(337);
hold on;
plot(time, L, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$L(t)$', 'Interpreter', 'latex');
title('Perdulet', 'Interpreter', 'latex');


subplot(338);
hold on;
plot(time, E_ossz, 'b', time, L, 'r');
xlabel('$time$', 'Interpreter', 'latex');
title('Osszenergia, perdulet', 'Interpreter', 'latex');
legend('Osszenergia', 'Perdulet','Location', 'east');


%% Masodrendu Runge-Kotta
%   

clc;
clear all;

%Allandok
m=5.972e+24;        %Fold tomege            (kg)
M=1.989e+30;        %Nap tomege             (kg)
G=6.67408e-11;      %Gravitacios allando    (m^3*kg^-1*s^-2)

%Kezdeti feltetelek
vx(1)=0;            %Napkozelben a sebesseg (m/s)
vy(1)=30.29e+3;     %Napkozelben a sebesseg (m/s)

rx(1)=147.1e+9;     %Napkozelben a tavolsag (m)
ry(1)=0;            %Napkozelben a tavolsag (m)

year=5;             %Hany evig keringjen a Fold a Nap korul
day=1;              %Lepeskoz, hany naponta szamoljon

t=60*60*24*365*year;%Keringesi ido          (secundum)
h=60*60*24*day;     %Lepeskoz               (secundum)
time(1)=0;          %Kezdettol eltelt ido   (secundum)

%Szamitas Masodrendu Runge-Kotta 

E_mozg(1)=(vx(1)^2+vy(1)^2)/2*m;            %Mozgasi energia
E_hely(1)=-G*M*m/(rx(1)^2+ry(1)^2)^0.5;     %Helyzeti energia
E_ossz(1)=E_mozg(1)+E_hely(1);              %Osszenergia

L(1)=m*norm(cross([rx(1),ry(1),0],[vx(1),vy(1),0]));    %Perdulet

for i = 2:t/h
    vxk1 = -G*M*rx(i-1)*h/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    vyk1 = -G*M*ry(i-1)*h/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    rxk1 = vx(i-1)*h;
    ryk1 = vy(i-1)*h;

    vxk2 = -G*M*h*(rx(i-1)+rxk1/2)/((rx(i-1)+rxk1/2)^2+(ry(i-1)+ryk1/2)^2)^(3/2);
    vyk2 = -G*M*h*(ry(i-1)+ryk1/2)/((rx(i-1)+rxk1/2)^2+(ry(i-1)+ryk1/2)^2)^(3/2);
    rxk2 = (vx(i-1)+vxk1/2)*h;
    ryk2 = (vy(i-1)+vyk1/2)*h;

    vx(i) = vx(i-1)+vxk2;
    vy(i) = vy(i-1)+vyk2;
    rx(i) = rx(i-1)+rxk2;
    ry(i) = ry(i-1)+ryk2;
    
    E_mozg(i)=(vx(i)^2+vy(i)^2)/2*m;
    E_hely(i)=-G*M*m/(rx(i)^2+ry(i)^2)^0.5;
    E_ossz(i)=E_mozg(i)+E_hely(i);
    
    L(i)=m*norm(cross([rx(i), ry(i), 0],[vx(i), vy(i), 0]));
    
    time(i)=time(i-1)+h;
end

%Kiiratas

figure;
str=sprintf('Masodrendu Runge-Kotta \n(keringesi ido %d ev, lepeskoz %d nap)', year, day);
sgtitle(str,'Interpreter', 'latex');
subplot(331);
hold on;
%axis equal;
plot(rx, ry, 'b');
plot(0, 0, 'oy', 'MarkerFaceColor', 'y');
xlabel('$r_x(t)$', 'Interpreter', 'latex');
ylabel('$r_y(t)$', 'Interpreter', 'latex');
title('Palyagorbe', 'Interpreter', 'latex');

subplot(332);
hold on;
plot(time, vx, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_x(t)$', 'Interpreter', 'latex');
title('X sebesseg-ido fv', 'Interpreter', 'latex');

subplot(333);
hold on;
plot(time, vy, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_y(t)$', 'Interpreter', 'latex');
title('Y sebesseg-ido fv', 'Interpreter', 'latex');

subplot(334);
hold on;
plot(time, E_mozg, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{mozg}(t)$', 'Interpreter', 'latex');
title('Mozgasi energia', 'Interpreter', 'latex');

subplot(335);
hold on;
plot(time, E_hely, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{hely}(t)$', 'Interpreter', 'latex');
title('Helyzeti energia', 'Interpreter', 'latex');

subplot(336);
hold on;
plot(time, E_ossz, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{ossz}(t)$', 'Interpreter', 'latex');
title('Osszenergia', 'Interpreter', 'latex');

subplot(337);
hold on;
plot(time, L, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$L(t)$', 'Interpreter', 'latex');
title('Perdulet', 'Interpreter', 'latex');


subplot(338);
hold on;
plot(time, E_ossz, 'b', time, L, 'r');
xlabel('$time$', 'Interpreter', 'latex');
title('Osszenergia, perdulet', 'Interpreter', 'latex');
legend('Osszenergia', 'Perdulet','Location', 'east');

%% Negyedrendu Runge-Kotta
%   

clc;
clear all;

%Allandok
m=5.972e+24;        %Fold tomege            (kg)
M=1.989e+30;        %Nap tomege             (kg)
G=6.67408e-11;      %Gravitacios allando    (m^3*kg^-1*s^-2)

%Kezdeti feltetelek
vx(1)=0;            %Napkozelben a sebesseg (m/s)
vy(1)=30.29e+3;     %Napkozelben a sebesseg (m/s)

rx(1)=147.1e+9;     %Napkozelben a tavolsag (m)
ry(1)=0;            %Napkozelben a tavolsag (m)

year=7;             %Hany evig keringjen a Fold a Nap korul
day=1;              %Lepeskoz, hany naponta szamoljon

t=60*60*24*365*year;%Keringesi ido          (secundum)
h=60*60*24*day;     %Lepeskoz               (secundum)
time(1)=0;          %Kezdettol eltelt ido   (secundum)

%Szamitas Negyedrendu Runge-Kotta 

E_mozg(1)=(vx(1)^2+vy(1)^2)/2*m;            %Mozgasi energia
E_hely(1)=-G*M*m/(rx(1)^2+ry(1)^2)^0.5;     %Helyzeti energia
E_ossz(1)=E_mozg(1)+E_hely(1);              %Osszenergia

L(1)=m*norm(cross([rx(1),ry(1),0],[vx(1),vy(1),0]));    %Perdulet

for i = 2:t/h
    vxk1 = -h*G*M*rx(i-1)/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    vyk1 = -h*G*M*ry(i-1)/(rx(i-1)^2+ry(i-1)^2)^(3/2);
    rxk1 = h*vx(i-1);
    ryk1 = h*vy(i-1);

    vxk2 = -h*G*M*(rx(i-1)+rxk1/2)/((rx(i-1)+rxk1/2)^2+(ry(i-1)+ryk1/2)^2)^(3/2);
    vyk2 = -h*G*M*(ry(i-1)+ryk1/2)/((rx(i-1)+rxk1/2)^2+(ry(i-1)+ryk1/2)^2)^(3/2);
    rxk2 = h*(vx(i-1)+vxk1/2);
    ryk2 = h*(vy(i-1)+vyk1/2);

    vxk3 = -h*G*M*(rx(i-1)+rxk2/2)/((rx(i-1)+rxk2/2)^2+(ry(i-1)+ryk2/2)^2)^(3/2);
    vyk3 = -h*G*M*(ry(i-1)+ryk2/2)/((rx(i-1)+rxk2/2)^2+(ry(i-1)+ryk2/2)^2)^(3/2);
    rxk3 = h*(vx(i-1)+vxk2/2);
    ryk3 = h*(vy(i-1)+vyk2/2);

    vxk4 = -h*G*M*(rx(i-1)+rxk3)/((rx(i-1)+rxk3)^2+(ry(i-1)+ryk3)^2)^(3/2);
    vyk4 = -h*G*M*(ry(i-1)+ryk3)/((rx(i-1)+rxk3)^2+(ry(i-1)+ryk3)^2)^(3/2);
    rxk4 = h*(vx(i-1)+vxk3);
    ryk4 = h*(vy(i-1)+vyk3);

    vx(i) = vx(i-1)+(vxk1+2*vxk2+2*vxk3+vxk4)/6;
    vy(i) = vy(i-1)+(vyk1+2*vyk2+2*vyk3+vyk4)/6;
    rx(i) = rx(i-1)+(rxk1+2*rxk2+2*rxk3+rxk4)/6;
    ry(i) = ry(i-1)+(ryk1+2*ryk2+2*ryk3+ryk4)/6;
    
    E_mozg(i)=(vx(i)^2+vy(i)^2)/2*m;
    E_hely(i)=-G*M*m/(rx(i)^2+ry(i)^2)^0.5;
    E_ossz(i)=E_mozg(i)+E_hely(i);
    
    L(i)=m*norm(cross([rx(i), ry(i), 0],[vx(i), vy(i), 0]));
    
    time(i)=time(i-1)+h;
end

%Kiiratas

figure;
str=sprintf('Negyedrendu Runge-Kotta \n(keringesi ido %d ev, lepeskoz %d nap)', year, day);
sgtitle(str,'Interpreter', 'latex');
subplot(331);
hold on;
%axis equal;
plot(rx, ry, 'b');
plot(0, 0, 'oy', 'MarkerFaceColor', 'y');
xlabel('$r_x(t)$', 'Interpreter', 'latex');
ylabel('$r_y(t)$', 'Interpreter', 'latex');
title('Palyagorbe', 'Interpreter', 'latex');

subplot(332);
hold on;
plot(time, vx, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_x(t)$', 'Interpreter', 'latex');
title('X sebesseg-ido fv', 'Interpreter', 'latex');

subplot(333);
hold on;
plot(time, vy, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_y(t)$', 'Interpreter', 'latex');
title('Y sebesseg-ido fv', 'Interpreter', 'latex');

subplot(334);
hold on;
plot(time, E_mozg, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{mozg}(t)$', 'Interpreter', 'latex');
title('Mozgasi energia', 'Interpreter', 'latex');

subplot(335);
hold on;
plot(time, E_hely, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{hely}(t)$', 'Interpreter', 'latex');
title('Helyzeti energia', 'Interpreter', 'latex');

subplot(336);
hold on;
plot(time, E_ossz, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{ossz}(t)$', 'Interpreter', 'latex');
title('Osszenergia', 'Interpreter', 'latex');

subplot(337);
hold on;
plot(time, L, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$L(t)$', 'Interpreter', 'latex');
title('Perdulet', 'Interpreter', 'latex');

subplot(338);
hold on;
plot(time, E_ossz, 'b', time, L, 'r');
xlabel('$time$', 'Interpreter', 'latex');
title('Osszenergia, perdulet', 'Interpreter', 'latex');
legend('Osszenergia', 'Perdulet','Location', 'east');

%%
%   Ode45

clc;
clear all;

%Allandok
m=5.972e+24;        %Fold tomege            (kg)
M=1.989e+30;        %Nap tomege             (kg)
G=6.67408e-11;      %Gravitacios allando    (m^3*kg^-1*s^-2)

%Kezdeti feltetelek
vx_0=0;            %Napkozelben a sebesseg (m/s)
vy_0=30.29e+3;     %Napkozelben a sebesseg (m/s)

rx_0=147.1e+9;     %Napkozelben a tavolsag (m)
ry_0=0;            %Napkozelben a tavolsag (m)

year=6;            %Hany evig keringjen a Fold a Nap korul

t=60*60*24*365*year;%Keringesi ido          (secundum)

%Szamitas Ode45

tspan = [0 t];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
kezdeti_feltetel=[vx_0; vy_0; rx_0; ry_0];
fv=@(time,r)    [                           % r(1)=vx, r(2)=vy, r(3)=rx, r(4)=ry
    -G*M*r(3)/(r(3)^2+r(4)^2)^1.5;
    -G*M*r(4)/(r(3)^2+r(4)^2)^1.5;
    r(1);
    r(2);
    ];
[time,r] = ode45(fv,tspan,kezdeti_feltetel,options);

E_mozg=m*( r(:,1).^2 + r(:,2).^2 )/2;       %Mozgasi energia
E_hely=-G*M*m/(r(:,3).^2+r(:,4).^2).^0.5;   %Helyzeti energia
E_ossz=E_mozg+E_hely;                       %Osszenergia

%L=m*norm(cross([r(:,3),r(:,4),zeros(length(r),1)],[r(:,1),r(:,2),zeros(length(r),1)]));    %Perdulet

for i=1:length(r)
    L(i)= m*norm(cross([r(i,3),r(i,4),0],[r(i,1),r(i,2),0]));
end    %Perdulet

%Kiiratas

figure;
str=sprintf('Ode45 \n(keringesi ido %d ev)', year);
sgtitle(str,'Interpreter', 'latex');
subplot(331);
hold on;
%axis equal;
plot(r(:,3), r(:,4), 'b');
plot(0, 0, 'oy', 'MarkerFaceColor', 'y');
xlim([-2e+11 2e+11]);
ylim([-2e+11 2e+11]);
xlabel('$r_x(t)$', 'Interpreter', 'latex');
ylabel('$r_y(t)$', 'Interpreter', 'latex');
title('Ode45', 'Interpreter', 'latex');

subplot(332);
hold on;
plot(time, r(:,1), 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_x(t)$', 'Interpreter', 'latex');
title('X sebesseg-ido fv', 'Interpreter', 'latex');

subplot(333);
hold on;
plot(time, r(:,2), 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$v_y(t)$', 'Interpreter', 'latex');
title('Y sebesseg-ido fv', 'Interpreter', 'latex');

subplot(334);
hold on;
plot(time, E_mozg, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{mozg}(t)$', 'Interpreter', 'latex');
title('Mozgasi energia', 'Interpreter', 'latex');

subplot(335);
hold on;
plot(time, E_hely, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{hely}(t)$', 'Interpreter', 'latex');
title('Helyzeti energia', 'Interpreter', 'latex');

subplot(336);
hold on;
plot(time, E_ossz, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$E_{ossz}(t)$', 'Interpreter', 'latex');
title('Osszenergia', 'Interpreter', 'latex');

subplot(337);
hold on;
plot(time, L, 'b');
xlabel('$time$', 'Interpreter', 'latex');
ylabel('$L(t)$', 'Interpreter', 'latex');
title('Perdulet', 'Interpreter', 'latex');

subplot(338);
hold on;
plot(time, E_ossz, 'b', time, L, 'r');
xlabel('$time$', 'Interpreter', 'latex');
title('Osszenergia, perdulet', 'Interpreter', 'latex');
legend('Osszenergia','Perdulet','Location', 'east');

%% Polarkoordinatak
%   

clc;
clear all;

%Allandok
m=5.972e+24;        %Fold tomege            (kg)
M=1.989e+30;        %Nap tomege             (kg)
G=6.67408e-11;      %Gravitacios allando    (m^3*kg^-1*s^-2)

%Kezdeti feltetelek
vx_0=0;             %Napkozelben a sebesseg (m/s)
vy_0=30.29e+3;      %Napkozelben a sebesseg (m/s)

rx_0=147.1e+9;      %Napkozelben a tavolsag (m)
ry_0=0;             %Napkozelben a tavolsag (m)

year=11;           %Hany evig keringjen a Fold a Nap korul

t=60*60*24*365*year;%Keringesi ido          (secundum)

%Szamitas Ode45 (polarkoordinatak)
%r..-h^2/r^3=-GM/r^2; h=r^2*fi., azaz a perdulet

tspan = [0 t];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
r_0 = (rx_0^2+ry_0^2)^(1/2);
h = r_0*(vx_0^2+vy_0^2)^0.5;                %r^2*teta derivalt, perdulet egysegnyi tomegre
kezdeti_feltetel=[r_0, 0, 0];               %r, kezdeti teta derivalt,  teta
fv=@(time,r)    [                           
    r(2);    
    -G*M/r(1)^2+h^2/r(1)^3;    
    h/r(1)^2;
    ];
[time,r] = ode45(fv,tspan,kezdeti_feltetel,options);

E_mozg=m*r(:,1).^2/2;                       %Mozgasi energia
E_hely=-G*M*m/r(:,1);                       %Helyzeti energia
E_ossz=E_mozg+E_hely;                       %Osszenergia

L= m.*r(:,1).*r(:,2);

figure;
hold on;
plot(cos(r(:,3)).*r(:, 1), sin(r(:,3)).*r(:,1));
plot(0, 0, 'oy', 'MarkerFaceColor', 'y');
xlim([-2e+11 2e+11]);
ylim([-2e+11 2e+11]);
xlabel('$r_x(t)$', 'Interpreter', 'latex');
ylabel('$r_y(t)$', 'Interpreter', 'latex');
str=sprintf('Ode45 polarkoordinatakkal \n(keringesi ido %d ev)', year);
title(str,'Interpreter', 'latex');