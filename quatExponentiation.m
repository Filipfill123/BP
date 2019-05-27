function Q_t = quatExponentiation(Q,t)
% Q_t = Q^t = exp(t * log(Q)), t in R

Q_t = quatExp(t*quatLog(Q));