function [alphan,ok]=Wolfe(alpha,x,D,Oracle)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1;
   omega2 = 0.9;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F0,G0] = Oracle(x,ind);

   // Initialisation de l'algorithme

   alphan = alpha;
   xn     = x;

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

   while ok == 0
       
      
      xp = xn;
      xn = x + (alphan*D);

      // Calcul des conditions de Wolfe

      // 1ere condition de Wolfe
      cond1=%F;
      cond2=%F;

      [Fn,Gn] = Oracle(xn,4);
      if (Fn <= (F0 + omega1*alphan*(G0'*D))) then
        cond1 = %T;
      end
      
      // 2ere condition de Wolfe
      if (Gn'*D) >= omega2*(G0'*D) then
        cond2 = %T;
      end      
      
      
      // -----> A completer...

      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      if cond1&cond2==%T then
        ok = 1
      end
      
      // - sinon, modifier la valeur de alphan : on reboucle.
 
      if cond1 == %F then
        alphamax = alphan;
        alphan = (1/2)*(alphamax + alphamin);
      else
          if cond2 == %F then
            alphamin = alphan;
            if isinf(alphamax) then
              alphan = (2)*(alphamin);
            else
              alphan = (1/2)*(alphamax + alphamin);
            end
          else
            ok=1;
          end
      end
      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end

endfunction
