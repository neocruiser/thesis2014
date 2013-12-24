#########################################################################
#								SECTION -1-								#
#########################################################################
info: http://tinyurl.com/pwv7q38

## GET PROFILE DATA ON A SCHOLAR
library(scholar)
tremblay <- 'qThJPrgAAAAJ&hl'
moraga <- 'TNm2HpwAAAAJ&hl'
bassim <- 'mDJZY3oAAAAJ'
profileUser1 <- get_profile(tremblay)
profileUser2 <- get_profile(moraga)
profileUser3 <- get_profile(bassim)

library(ggplot2)
ids <- c(tremblay, moraga)
compareUsers <- compare_scholar_careers(ids)
ggplot(compareUsers, aes(x=career_year, y=cites)) +
geom_line(aes(linetype=name)) + theme_bw()


