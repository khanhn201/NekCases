% load("ha10n10.mat")
format long;

[z,~]=zwgll(N);

nelxy = length(meshVel) - 1;
nelz = 1;

z_length = 10;

delta_el = diff(meshVel);
el_center = (meshVel(2:end) + meshVel(1:end-1))/2.0;

n2 = (N+1)^2;
n3 = (N+1)^3;

data_r    = zeros(nelxy*nelxy*nelz, n3, 6); % xyzuvw
data_i    = zeros(nelxy*nelxy*nelz, n3, 6); % xyzuvw
data_base = zeros(nelxy*nelxy*nelz, n3, 6); % xyzuvw
for j=1:nelxy
for i=1:nelxy
    el = i + (j-1)*nelxy;
    idxV  = (el-1)*n2 + (1:n2);

    zfx =  el_center(i) + z*delta_el(i)/2.0;
    zfy =  el_center(j) + z*delta_el(j)/2.0;
    zfz =  z_length/2 + z*z_length/2.0; % [0, z_length]

    uel = reshape(  uvec(idxV), N+1, N+1);
    vel = reshape(  vvec(idxV), N+1, N+1);
    wel = reshape(  wvec(idxV), N+1, N+1);
    Uel = reshape(U_base(idxV), N+1, N+1); % (Nf,Nf)
    uel = repmat(uel, 1, 1, N+1);
    vel = repmat(vel, 1, 1, N+1);
    wel = repmat(wel, 1, 1, N+1);
    Uel = repmat(Uel, 1, 1, N+1);     % (Nf,Nf,Nf)
    uel = uel(:);
    vel = vel(:);
    wel = wel(:);
    Uel = Uel(:); % (Nf^3, 1)



    [X,Y,Z] = ndgrid(zfx,zfy,zfz);

    data_r(el,:,1) = X(:);
    data_r(el,:,2) = Y(:);
    data_r(el,:,3) = Z(:);
    data_i(el,:,1) = X(:);
    data_i(el,:,2) = Y(:);
    data_i(el,:,3) = Z(:);
    data_base(el,:,1) = X(:);
    data_base(el,:,2) = Y(:);
    data_base(el,:,3) = Z(:);
    data_r(el,:,4) = real(uel);
    data_r(el,:,5) = real(vel);
    data_r(el,:,6) = real(wel);
    data_i(el,:,4) = imag(uel);
    data_i(el,:,5) = imag(vel);
    data_i(el,:,6) = imag(wel);
    data_base(el,:,6) = Uel;
end
end





E = nelxy*nelxy*nelz;
lr1 = [N+1,N+1,N+1];
elmap = (1:E);
istep = 0;
time= 0;
fields = ['XU'];
emode = 'le';
wdsz = 8; % precision
etag = 6.54321;

status = writenek('u_real0.f00001',
                  data_r,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
status = writenek('u_imag0.f00001',
                  data_i,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
status = writenek('u_base0.f00001',
                  data_base,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
