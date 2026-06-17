function [Q,gid] = set_tp_semq(Nelx,Nely,N)

N1 = N+1;

Nxg = Nelx*N + 1;
Nyg = Nely*N + 1;

Ng = Nxg * Nyg;

E = Nelx*Nely;

gid = zeros(N1,N1,E);

e = 0;

for ey = 1:Nely
for ex = 1:Nelx

    e = e + 1;

    gx = (ex-1)*N + (1:N1);
    gy = (ey-1)*N + (1:N1);

    for j=1:N1
    for i=1:N1

        gid(i,j,e) = gx(i) + (gy(j)-1)*Nxg;

    end
    end

end
end

Q = sparse(E*N1*N1, Ng);

row = 0;

for e=1:E
for j=1:N1
for i=1:N1

    row = row + 1;

    col = gid(i,j,e);

    Q(row,col) = 1;

end
end
end

