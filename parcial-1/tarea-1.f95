program Punto1
     integer :: i
     real, parameter :: pi=3.14159
     real, parameter :: e=2.71828
     double precision, parameter :: h=6.62607e-34
     double precision, parameter :: kb=1.38065e-23
     double precision :: T, V, G, U, S, N
     double precision, dimension(8) :: m=(/6.64e-27, 3.35e-26, 4.49e-26, 2.28e-25, 6.63e-26, 1.49e-26, 1.32e-25, 1.99e-26/) !Vector para masa atomica de cada elemento de la tabla
     double precision, dimension(8) :: Qe=(/1.00, 1.00, 2.00, 1.00, 1.00, 1.00, 2.00, 1.00/) !Vector para estado de ge1 de cada elemento de la tabla
     double precision:: Gt(8), Ut(8), St(8), Ge(8), Ue(8), Se(8), Gte(8), Ute(8), Ste(8)
     !m(1)=Helio m(2)=Neon m(3)=Aluminio m(4)=Bario m(5)=Argon m(6)=Berilio m(7)=Bromo m(8)=Carbono
     T=298
     N=6.023e23
     V=4.11e-26 !Volumen calculado por ley de gases ideales y en m^3
   
     !Calculo con la funcion de particion traslacional
     open(2,file = 'Funcion_Traslacional.dat',status="unknown")
     write(2,*) "          G", "                          U", "                         S"
     Do i=1, 8
       G=8.43e-5*T*(-log((((2*pi*m(i)*kb*T)/(h**2))**(1.5))*V)+1)!Energia libre de Gibs
       U=N*kb*T*1.5/1000 !Energia interna
       S=N*((1.5*kb)+(kb*log((((2*pi*m(i)*kb*T)/(h**2))**(1.5))*V))) !Entropia
       write(2,*) G, U, S
       Gt(i)=G
       Ut(i)=u
       St(i)=S
       G=0
       U=0
       S=0
     end Do  
     close(3)
     G=0
     U=0
     S=0
     !Calculo con la funcion de particion electronica
     open(3,file = 'Funcion_Electronica.dat',status="unknown")
     write(3,*) "          G", "                          U", "                         S"
     Do i=1, 8
       G=((-kb*T*log(Qe(i)))+(kb*T)) !Energia libre de Gibs
       U=kb*(T**2)*(1/Qe(i))*10e9 !Energia interna
       S=N/1000*kb*T*((1/Qe(i))+(kb*log(Qe(i)))) !Entropia
       write(3,*) G, U, S
       Ge(i)=G
       Ue(i)=u
       Se(i)=S
       G=0
       U=0
       S=0
     end Do 
     close(2)
   !Calculo con la funcion de particion electronica y traslacional
     open(4,file = 'Funcion_Traslacional+Electronica.dat',status="unknown")
     write(4,*) "          G", "                          U", "                         S"
     Do i=1, 8
       G=8.31e-3*T*log((((((2*pi*m(i)*kb*T)/(h**2))**(1.5))*V)*Qe(i))/N)!Energia libre de Gibs
       U=Ue(i)+Ut(i) !Energia interna
       S=Se(i)+St(i)!Entropia
       write(4,*) G, U, S
       Gte(i)=G
       Ute(i)=u
       Ste(i)=S
       G=0
       U=0
       S=0
     end Do 
     close(4)
     end program Punto1