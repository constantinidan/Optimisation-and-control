function [fopt,xopt,gopt]=Gradient_PR(OraclePG,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient Polak-Ribiere                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres du gradient a pas fixe";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"1";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   [F,G] = OraclePG(x,4);
   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient
      Fp = F
      Gp = G
      [F,G] = OraclePG(x,4);

//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de beta
      if k == 1 then
        D = -G;
      else
        Beta = (G')*(G-Gp)/(Gp'*Gp)
    
//    - calcul de la direction de descente
        D = -G + Beta*Dp;
      end
      Dp = D
      
//    - calcul de la longueur du pas de gradient

      [alpha,ok]=Wolfe(alphai,x,D,OraclePG)

//    - mise a jour des variables

      x = x + (alpha*D);

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de gradient a pas fixe')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction



