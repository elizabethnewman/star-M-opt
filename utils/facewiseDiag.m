function[d] = facewiseDiag(A,k)
    % starM product for third-order tensors 
    %
    % Inputs:
    %   A  : n1 x n2 x n3 array or A : r x 1 x n3 where r is the number of entries on the kth diagonal
    %   k : the kth diagonal (optional). k=0 represents the main diagonal (default), k>0 is above the main diagonal, and k<0 is below the main diagonal.
    %
    % Outputs:
    %   d : r x n3 or d : n1 x n2 x n3
    %

    if ~exist('k','var') || isempty(k), k = 0; end

    if size(A,1) == 1 || size(A,2) == 1
        % outputs tensor with square frontal slices
        % input is numel(idx) x 1 x tubesize or 1 x numel(idx) x tubesize

        % output the diagonal tubes of the tensor
        n   = max(size(A,1),size(A,2));
        szA = size(A);
        nd  = prod(szA(3:end));


        % main diagonal
        idx0 = (0:(n - 1)) * n + (1:n);
    
        % change which diagonal
        idx = idx0 - k;
        idx = idx(idx >= 1 & idx <= n^2);

        % adjust for multiple frontal slices
        idx2 = (0:(nd - 1)) * n^2;

        % form main index
        IDX = kron(ones(1,nd),idx) + kron(idx2,ones(1,numel(idx)));

        % output tensor
        d = zeros([n,n,szA(3:end)]);
        % form tensor
        d(IDX) = A(:)';
    else
        % output the diagonal tubes of the tensor
        m  = size(A,1);
        n  = size(A,2);
        szA = size(A);
        nd = prod(szA(3:end));
        
        r = min(m,n);
    
        % main diagonal
        idx0 = (0:(r - 1)) * m + (1:r);
    
        % change which diagonal
        idx = idx0 - k;
        idx = idx(idx >= 1 & idx <= m * n);
        
        % adjust for multiple frontal slices
        idx2 = (0:(nd - 1)) * (m * n);
    
        % form main index
        IDX = kron(ones(1,nd),idx) + kron(idx2,ones(1,numel(idx)));
    
        % output diag
        d = reshape(A(IDX),[numel(idx),1,szA(3:end)]);
    end
end

