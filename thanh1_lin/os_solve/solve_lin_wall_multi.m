function [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall_multi(N, meshVel, meshWall, rho, mu, sigma, sigma_w, W, U,Phi, By, alpha)
  meshPhi  = [-meshWall(end:-1:1), meshVel(2:end-1), meshWall];

  NelWall  = length(meshWall) - 1;
  NelxV    = length(meshVel)  - 1;
  NelyV    = length(meshVel)  - 1;
  NelxP    = length(meshPhi)  - 1;
  NelyP    = length(meshPhi)  - 1;

  NelV     = NelxV*NelyV;
  NelP     = NelxP*NelyP;

  [Ah,Bh,~,Dh,z,~] = semhat(N);
  n = N + 1;
  n2 = n^2;
  Ih = speye(n); 

  [QVel, ~] = set_tp_semq(NelxV,NelyV,N);
  [QPhi, ~] = set_tp_semq(NelxP,NelyP,N);
  [QPr , ~] = set_tp_semq(NelxV,NelyV,N-2);
  QPr = speye(size(QPr,1));

  IhGlo = speye(n*NelyV - (NelyV-1));
  R = IhGlo(2:end-1,:);
  RVel = kron(R, R);           % zero dirichlet for velocity
  RPhi = speye(size(QPhi,2));  % zero neumann   for phi

  [zd,wd]=zwgll(ceil(3*N/2)); Bhd=diag(wd); % For dealiasing
  JD    = interp_mat(zd,z);
  JDXY  = kron(JD ,JD );
  Ihd   = speye(size(zd,1));
  Dhd   = deriv_mat(zd);



  [zp,wp]=zwgl(N-1); Bhp=diag(wp); % For pressure
  np2 = (n-2)^2;
  JP    = interp_mat(zp,z);
  JPT   = interp_mat(z,zp);
  JPXY  = kron(JP ,JP );
  JPTXY = kron(JPT,JPT);
  Ihp   = speye(size(zp,1));


  Kuu     = zeros(n2*NelV,n2*NelV);
  Kvv     = zeros(n2*NelV,n2*NelV);
  Kwu     = zeros(n2*NelV,n2*NelV);
  Kwv     = zeros(n2*NelV,n2*NelV);
  Kww     = zeros(n2*NelV,n2*NelV);
  Kphiphi = zeros(n2*NelP,n2*NelP);
  Kphiu   = zeros(n2*NelP,n2*NelV);
  Kphiv   = zeros(n2*NelP,n2*NelV);
  Kphiw   = zeros(n2*NelP,n2*NelV);
  Kuphi   = zeros(n2*NelV,n2*NelP);
  Kvphi   = zeros(n2*NelV,n2*NelP);
  Kwphi   = zeros(n2*NelV,n2*NelP);

  Kup     = zeros(n2*NelV,np2*NelV);
  Kvp     = zeros(n2*NelV,np2*NelV);
  Kwp     = zeros(n2*NelV,np2*NelV);
  Kpu     = zeros(np2*NelV,n2*NelV);
  Kpv     = zeros(np2*NelV,n2*NelV);
  Kpw     = zeros(np2*NelV,n2*NelV);

  Mu       = zeros(n2*NelV,1);
  Mv       = zeros(n2*NelV,1);
  Mw       = zeros(n2*NelV,1);

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
      
      Ahx  = Ah /(hx/2); Ahy  = Ah /(hy/2);
      Bhx  = Bh *(hx/2); Bhy  = Bh *(hy/2);
      Bhpx = Bhp*(hx/2); Bhpy = Bhp*(hy/2);
      Bhdx = Bhd*(hx/2); Bhdy = Bhd*(hy/2);
      
      AhXY  = kron(Bhy,Ahx) + kron(Ahy,Bhx);
      BhXY  = kron(Bhy,Bhx);
      BhpXY = kron(Bhpy,Bhpx);
      BhdXY = kron(Bhdy,Bhdx);

      DhX  = kron(Ih ,Dh )/(hx/2);
      DhY  = kron(Dh ,Ih )/(hy/2);
      DhdX = kron(Ihd,Dhd)/(hx/2);
      DhdY = kron(Dhd,Ihd)/(hy/2);

      eP   = (ey-1)*NelxP + ex;
      idxP = (eP-1)*n2 + (1:n2);
      if inFluid
          exV = ex - NelWall;
          eyV = ey - NelWall;

          eV    = (eyV-1)*NelxV + exV;
          idxV  = (eV-1)*n2 + (1:n2);
          idxPr = (eV-1)*np2 + (1:np2);

          Uel = JDXY*U(idxV); % Interp mean flow to dealias
          Kuu(idxV,idxV) = mu*AhXY + mu*alpha^2*BhXY...
                           + 1i*alpha*rho*JDXY'*BhdXY*diag(     Uel)*JDXY... % dealias
                           + sigma*By^2*BhXY;
          Kvv(idxV,idxV) = mu*AhXY + mu*alpha^2*BhXY...
                           + 1i*alpha*rho*JDXY'*BhdXY*diag(     Uel)*JDXY;   % dealias
          Kww(idxV,idxV) = mu*AhXY + mu*alpha^2*BhXY...
                           + 1i*alpha*rho*JDXY'*BhdXY*diag(     Uel)*JDXY... % dealias
                           + sigma*By^2*BhXY;

          Kwu(idxV,idxV) =            rho*JDXY'*BhdXY*diag(DhdX*Uel)*JDXY;
          Kwv(idxV,idxV) =            rho*JDXY'*BhdXY*diag(DhdY*Uel)*JDXY;


          Kuphi(idxV,idxP) = -1i*alpha*sigma*By     *BhXY;
          Kvphi(idxV,idxP) =                       0*BhXY;
          Kwphi(idxV,idxP) =          -sigma*By*DhX'*BhXY; % Weak form

          Kphiu(idxP,idxV) =  1i*alpha*sigma*By     *BhXY;
          Kphiv(idxP,idxV) =                       0*BhXY;
          Kphiw(idxP,idxV) =          -sigma*By     *BhXY*DhX;


          Kup(idxV,idxPr) =     -DhX'*BhXY*JPTXY; % Weak form
          Kvp(idxV,idxPr) =     -DhY'*BhXY*JPTXY; % Weak form
          Kwp(idxV,idxPr) =  1i*alpha*BhXY*JPTXY;

          Kpu(idxPr,idxV) =          -BhpXY*JPXY*DhX;
          Kpv(idxPr,idxV) =          -BhpXY*JPXY*DhY;
          Kpw(idxPr,idxV) = -1i*alpha*BhpXY*JPXY;

          Mu(idxV,1) = -rho*diag(BhXY);
          Mv(idxV,1) = -rho*diag(BhXY);
          Mw(idxV,1) = -rho*diag(BhXY);
      end
      Kphiphi(idxP,idxP) = sigmal*AhXY + sigmal*alpha^2*BhXY;
  end
  end




  KVelVel = [
        RVel*QVel'*Kuu*QVel*RVel', 0*RVel*QVel'*Kuu*QVel*RVel', 0*RVel*QVel'*Kuu*QVel*RVel';
      0*RVel*QVel'*Kvv*QVel*RVel',   RVel*QVel'*Kvv*QVel*RVel', 0*RVel*QVel'*Kvv*QVel*RVel';
        RVel*QVel'*Kwu*QVel*RVel',   RVel*QVel'*Kwv*QVel*RVel',   RVel*QVel'*Kww*QVel*RVel';
  ];
  KVelPhi = [
      RVel*QVel'*Kuphi*QPhi*RPhi';
      RVel*QVel'*Kvphi*QPhi*RPhi';
      RVel*QVel'*Kwphi*QPhi*RPhi';
  ];
  
  KVelPr = [
      RVel*QVel'*Kup*QPr;  
      RVel*QVel'*Kvp*QPr;  
      RVel*QVel'*Kwp*QPr;  
  ];
  KPrVel = [QPr'*Kpu*QVel*RVel', QPr'*Kpv*QVel*RVel', QPr'*Kpw*QVel*RVel'];
  KPhiVel = [RPhi*QPhi'*Kphiu*QVel*RVel', RPhi*QPhi'*Kphiv*QVel*RVel', RPhi*QPhi'*Kphiw*QVel*RVel'];
  KPhiPhi = RPhi*QPhi'*Kphiphi*QPhi*RPhi';

  Du = spdiags(Mu,0,numel(Mu),numel(Mu));
  Dv = spdiags(Mv,0,numel(Mv),numel(Mv));
  Dw = spdiags(Mw,0,numel(Mw),numel(Mw));
  MV = blkdiag(RVel*QVel'*Du*QVel*RVel',RVel*QVel'*Dv*QVel*RVel',RVel*QVel'*Dw*QVel*RVel');




  KPhiPhi = sparse(KPhiPhi);
  KPhiVel = sparse(KPhiVel);

  disp('invert')
  sol = KPhiPhi\KPhiVel;
  disp('invert2')
  KVV2 = KVelVel - KVelPhi*sol;
  [KL,KU,KP] = lu(KVV2);
  disp('lu')
  KVP2 = KU \ (KL \ (KP*KVelPr));
  MV2 = KU \ (KL \ (KP*MV));
  M2 = (KPrVel*KVP2) \ (KPrVel*MV2);
  disp('done invert')

  M3 = MV - KVelPr*M2;
  % M4 = KU \ (KL \ (KP*M3));
  M4 = MV2 - KVP2*M2;

  disp('eigen')
  [vecs, gamma] = eig(M4);
  % [vecs, gamma] = eigs(M3, KVV2, 1000);
  gamma = diag(gamma);
  gamma = 1./gamma;
  disp('done eigen')


  KV = KVV2 * vecs;
  MV = M3   * vecs;
  R = KV - MV * diag(gamma);
  res = vecnorm(R) ./ (vecnorm(KV) + vecnorm(MV) + eps);

  good = res < 1e-8;

  gamma = gamma(good);
  vecs  = vecs(:,good);
  res   = res(good);
  [re_sorted, idx] = sort(real(gamma), 'descend');
  disp(gamma(idx(1:min(5,end))))

  npos = nnz(real(gamma) > 0);
  if npos == 0
      warning('No unstable eigenvalues found (Re(gamma) > 0).');
  elseif npos > 1
      warning('More than one unstable eigenvalues found (Re(gamma) > 0).');
  end

  [~, unstable] = max(real(gamma));
  gamma = gamma(unstable);
  vu = vecs(:, unstable);
  vdof = size(RVel, 1);
  uvec   =  QVel*RVel'*vu(0*vdof+1:1*vdof);
  vvec   =  QVel*RVel'*vu(1*vdof+1:2*vdof);
  wvec   =  QVel*RVel'*vu(2*vdof+1:3*vdof);
  pvec   =  QPr*gamma*M2*vu;
  phivec = -QPhi*RPhi'*sol*vu;
end
