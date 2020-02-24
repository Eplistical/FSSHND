
function hamiltonian()
    w = 1.0;
    ww = 2.0;
    x0 = 1.0;

    m = 1.0;
    Vc = 0.2;
    W = 1.0;
    N = 2;

    function Vk = cal_V(x, k)
        Vk = 0.5 * m * w^2 * (sum(x.^2) - x(k)^2) + 0.5 * m * ww^2 * (x(k) - x0)^2;
    end

    function H = cal_H(x)
        H = zeros(N,N);
        % site energy
        for k = 1:N
            H(k,k) = cal_V(x, k);
        end
        % neighbor sites coupling
        for k = 1:N-1
            H(k, k+1) = Vc * exp(1i * W * (x(k) + x(k+1)));
            H(k+1, k) = conj(H(k, k+1));
        end
    end

    Htest = cal_H([1,2])

    Nx = 50;
    xarr = linspace(-3,3,Nx);
    yarr = linspace(-3,3,Nx);
    [meshx, meshy] = meshgrid(xarr, yarr);
    E1 = zeros(Nx, Nx);
    E2 = zeros(Nx, Nx);
    for ix=1:Nx
        for iy=1:Nx
            x = meshx(ix, iy);
            y = meshy(ix, iy);
            h = cal_H([x;y]);
            [evt, eva] = eig(h);
            E1(ix,iy) = eva(1,1);
            E2(ix,iy) = eva(2,2);
        end
    end
    contour(meshx, meshy, E1, 50);
end
