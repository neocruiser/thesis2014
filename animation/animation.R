require(animation)

ani.options(interval = 0.3)
lln.ani(FUN = function(n, mu) rchisq(n, df = mu), mu = 5, cex = 0.6)

## Lilac chaser
ani.options(nmax = 20)
par(mar = c(1, 1, 1, 1))
vi.lilac.chaser()