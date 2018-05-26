function [cm,cl,cp]=DiscreteVortexMethod(m,p,f,N,c,uinf,alfa,xref)
m=m*0.01; %combadura máxima
p=p*0.1; %x posicion combadura maxima
f=f*0.01; %espesor
i_p = 1;%posición matricial de la combadura máxima
min=c; %variable auxiliar utilizada para calcular la posición de i_p
% xref es el punto respecto el cual se quiere calcular el coeficiente de momentos
rho = 1; %De acuerdo a la nota alcarativa definimos la densidad unitaria, ya que no interviene en los coeficientes adimensionales.
sumatorio1=0;
for i = 1:N+1 %Hacemos el bucle para discretizar la cuerda del perfil
    x(i)=(c/2)*(1-cos((i-1)*(pi/(N+1))));
end

for i=1:N+1 %Buscamos que punto de nuestra discretizacion se acerca mas al de la combadura máxima
    if abs(x(i)-p*c)<min
        min=abs(x(i)-p*c);
        i_p=i;
    end
end

for i = 1:i_p %definimos la combadura del perfil hasta el máximo
    y(i)= m*(x(i)/(p^2))*(2*p-(x(i)/c));
end
for i = i_p:N+1 %la seguimos definiendo hasta el final
    y(i) = m*((c-x(i))/((1-p)^2))*(1+(x(i)/c)-2*p);
end

for i = 1:N %calculamos los deltas de x y de y
    delta_x(i)= x(i+1)-x(i);
    delta_y(i) = y(i+1)-y(i);
end

for i = 1:N %calculamos la longitud de los paneles y sus vectores normales y tangentes
    l_p(i) = sqrt((delta_x(i)^2)+(delta_y(i)^2)); %longitud del panel
    nx(i) = -delta_y(i)/l_p(i); %vector normal x. 
    ny(i) = delta_x(i)/l_p(i); %vector normal y
    tx(i) = delta_x(i)/l_p(i); %vector tangente x
    ty(i) = delta_y(i)/l_p(i); %vector tangente y
end

for i = 1:N %Buscamos los puntos importantes del panel
    vortex_pointx(i) = x(i) + delta_x(i)/4;
    vortex_pointy(i) = y(i) + delta_y(i)/4;
    control_pointx(i) = x(i) + (3/4)*delta_x(i);
    control_pointy(i) = y(i) + (3/4)*delta_y(i);
end

for i = 1:N %Calculamos la influencia de todos los vortices sobre cada uno de los puntos de control, asi como los coeficientes de influencia
    for j = 1:N
        radio2(i,j) = (control_pointx(i)-vortex_pointx(j))^2 + (control_pointy(i)-vortex_pointy(j))^2;
        u(i,j) = (1/(2*pi))*((control_pointy(i)-vortex_pointy(j))/(radio2(i,j))); %componente x de las velocidades indicidas
        w(i,j) = (-1/(2*pi))*((control_pointx(i)-vortex_pointx(j))/(radio2(i,j))); %componente y de las velocidades inducidas
        A(i,j) = u(i,j)*nx(i)+w(i,j)*ny(i); %coeficientes de influencia. Matriz A
    end
    RHS(i) = -uinf*(cos(alfa)*nx(i)+sin(alfa)*ny(i));  
end
RHS = RHS';
Gamma = inv(A)*RHS;

Cl = (2/uinf*c)*sum(Gamma); %Coeficiene de sustentación del perfil

for j = 1:N
    sumatorio = Gamma(j)*(vortex_pointx(j)-xref)*cos(alfa);
    sumatorio1=sumatorio+sumatorio1;
end
Cmref = ((-2)/(uinf*(c^2)))*sumatorio1; %Coeficiente de momentos respecto del punto de referencia seleccionado
for j = 1:N
    delta_Cp(j) = (2/uinf)*(Gamma(j)/l_p(j)); %Delta del coeficiente de momentos
end
for i=1:N
    z(i)=x(i);
end  
