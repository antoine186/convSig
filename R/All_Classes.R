setClass (
  # Class name
  "feat_obj",
  
  # Defining slot type
  representation (
    feat = "matrix",
    Multi_F = "matrix"
  )
)

setClass (
  # Class name
  "relu_obj",
  
  # Defining slot type
  representation (
    feat = "matrix",
    mat = "array",
    P = "matrix",
    LOSS = "numeric",
    test_LOSS = "numeric"
    
  )
)

setClass (
  # Class name
  "exp_obj",
  
  # Defining slot type
  representation (
    feat = "array",
    mat = "array",
    P = "matrix",
    LOSS = "numeric",
    test_LOSS = "numeric"
  )
)

setClass (
  # Class name
  "Shallowres",
  
  # Defining slot type
  representation (
    mut_mat = "matrix",
    wt = "numeric"
  )
)