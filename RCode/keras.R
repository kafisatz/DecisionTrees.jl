library(keras)

model <- keras_model_sequential() 
model %>% 
  layer_dense(units = q1, activation = 'tanh', input_shape = c(nrow(Xlearn))) %>% 
  layer_dense(units = 1, activation = k_exp)  

summary(model)

model %>% compile(
  loss = 'poisson',
  optimizer ='sgd',
)

fit <- model %>% fit(Xlearn, learn$ClaimNb, epochs=100, batch_size=10000)

