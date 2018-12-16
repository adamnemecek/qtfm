function Q = ddqft(q, mu, L)
% Discrete dual quaternion Fourier transform.
% To use this function it is necessary to install Quaternion and octonion
% toolbox for Matlab, see http://qtfm.sourceforge.net/

% Example: ddqft(randq(10,20,2), [quaternion(1,0,0) quaternion(0,1,0)], 'L')

% q - matrix of dual quaternions, in 1D case it has size (m, 1, 2), in 2D
% (m, n, 2);  q(:,:,1) is the real part of dual quaternions, q(:,:,2) is the dual part

% mu - dual quaternion, for which we have mu^2 = -1, so mu(1) must be pure
% unit quaternion and mu(2) is a quaternion orthogonal to mu(1)

% L = 'L' - left transform
% L = 'R' - right transform

Q = q;

f = 0.5*q(:,:,1);

if L == 'L'  % left transform
    f = mu(1)*mu(2)*f;
end

y = qfft2(f, -mu(1), L) + qfft2(-f, mu(1), L);

if L == 'R'  % right transform
    y = y*mu(1)*mu(2);
end

Q(:,:,1) = qfft2(q(:,:,1), mu(1), L);
Q(:,:,2) = qfft2(q(:,:,2), mu(1), L) + y;

end