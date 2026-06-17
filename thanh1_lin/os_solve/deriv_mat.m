      function[d] =  deriv_mat(x)

%     Compute the polynomial interpolation matrix from x to x

%     A bit faster and more stable for large n,  9/23/23 pff

      
      n = length(x);
      a = ones(n,1);

      s = 4/(max(x)-min(x)); % Allow large n

      e = ones(n,1);
      d = x*e' - e*x';

      for i=1:n;              %% Still has O(n^2) loop :/
        for j=1:(i-1); a(i)=s*a(i)*d(i,j); end;
        for j=(i+1):n; a(i)=s*a(i)*d(i,j); end;
      end;
      a=1./a; % These are the alpha_i's

      for j=1:n; d(j,j)=1; end;
      d=1./d;

      dd=zeros(n,1);
      for i=1:n; d(i,i)=0; dd(i)=sum(d(i,:)); end; % Row-sum = 0

      for j=1:n;  xj=x(j); aj=a(j);
         adx=a.*(x-xj); adx(j)=1;
         d(:,j) = aj./adx;
         d(j,j) = dd(j);
      end;

