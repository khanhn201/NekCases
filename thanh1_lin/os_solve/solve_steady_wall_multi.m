function [uvec,phivec,f] = solve_steady_wall_multi(N, meshVel, meshWall, mu, sigma, sigma_w, W, By)
  meshPhi  = [-meshWall(end:-1:1), meshVel(2:end-1), meshWall];

  NelWall  = length(meshWall) - 1;
  NelxV    = length(meshVel)  - 1;
  NelyV    = length(meshVel)  - 1;
  NelxP    = length(meshPhi)  - 1;
  NelyP    = length(meshPhi)  - 1;

  NelV     = NelxV*NelyV;
  NelP     = NelxP*NelyP;

  [Ah,Bh,~,Dh,~,~] = semhat(N);
  n = N + 1;
  n2 = n^2;
  Ih = speye(n); 

  [QVel, ~] = set_tp_semq(NelxV,NelyV,N);
  [QPhi, ~] = set_tp_semq(NelxP,NelyP,N);

  IhGlo = speye(n*NelyV - (NelyV-1));
  R = IhGlo(2:end-1,:);
  RVel = kron(R, R);           % zero dirichlet for velocity
  RPhi = speye(size(QPhi,2));  % zero neumann   for phi


  Kuu     = zeros(n2*NelV,n2*NelV);
  Kphiphi = zeros(n2*NelP,n2*NelP);
  Kphiu   = zeros(n2*NelP,n2*NelV);
  Kuphi   = zeros(n2*NelV,n2*NelP);
  rhs     = zeros(n2*NelV,1);

  for ey = 1:NelyP
  for ex = 1:NelxP
      inFluid = ...
          ex > NelWall && ex <= NelxP-NelWall && ...
          ey > NelWall && ey <= NelyP-NelWall;

      sigmal = sigma;
      if ~inFluid
         sigmal = sigma_w;
      end

      hx = meshPhi(ex+1) - meshPhi(ex);
      hy = meshPhi(ey+1) - meshPhi(ey);
      
      Bhx = Bh*(hx/2); Bhy = Bh*(hy/2);
      Ahx = Ah/(hx/2); Ahy = Ah/(hy/2);
      
      AhXY  = kron(Bhy,Ahx) + kron(Ahy,Bhx);
      BhXY  = kron(Bhy,Bhx);
      DhX   = kron(Ih ,Dh )/(hx/2);
      DhY   = kron(Dh ,Ih )/(hy/2);


      eP   = (ey-1)*NelxP + ex;
      idxP = (eP-1)*n2 + (1:n2);
      if inFluid
          exV = ex - NelWall;
          eyV = ey - NelWall;

          eV   = (eyV-1)*NelxV + exV;
          idxV = (eV-1)*n2 + (1:n2);

          Kuu  (idxV,idxV) = mu*AhXY + sigma*By^2*BhXY;
          Kuphi(idxV,idxP) = -sigma*By*DhX'*BhXY; % Weak form

          Kphiu(idxP,idxV) = -sigma*By*BhXY*DhX;
          rhs  (idxV,1)    = BhXY*ones(n2  , 1);
      end
      Kphiphi(idxP,idxP) = sigmal*AhXY;
  end
  end

  K = [
    RVel*QVel'*Kuu  *QVel*RVel', RVel*QVel'*Kuphi  *QPhi*RPhi';
    RPhi*QPhi'*Kphiu*QVel*RVel', RPhi*QPhi'*Kphiphi*QPhi*RPhi';
  ];
  
  rhsAll = [
    RVel*QVel'*rhs;
    RPhi*QPhi'*zeros(n2*NelP, 1);
  ];

  K = sparse(K);

  sol = K\rhsAll;
  % L = ichol(K);
  % [sol, flag, relres, iter] = pcg(K, rhsAll, 1e-8, 1000, L, L');
  % setup.type = 'ilutp';
  % setup.droptol = 1e-6;
  % [L,U] = ilu(K,setup);
  % restart = 50;
  % tol     = 1e-8;
  % maxit   = 2000;
  % [sol,flag,relres,iter,resvec] = gmres( ...
  %     K, rhsAll, restart, tol, maxit, L, U);

  % iter
  % relres

  uvec   = QVel*RVel'*sol(1:size(RVel,1));
  phivec = QPhi*RPhi'*sol(size(RVel,1)+1:end);

  umax = max(uvec);
  uvec = uvec/umax;
  phivec = phivec/umax;
  f = 1.0/umax;
end
