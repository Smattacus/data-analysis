function s = norm_gauss(x, y)
n = length( x );
x = reshape( x, n, 1 );
y = reshape( y, n, 1 );

%sort according to x
X = [x,y];
X = sortrows( X );
x = X(:,1);
y = X(:,2);

dx = diff( x );
dy = 0.5*(y(1:length(y)-1) + y(2:length(y)));
s = sum( dx .* dy );
if( s > 1.5 | s < 0.5 )
    fprintf( 'Data is not normalized! The pdf sums to: %f.\n\r', s );
end
