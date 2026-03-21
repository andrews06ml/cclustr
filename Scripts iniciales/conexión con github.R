getwd()


# Crear mi usuario
usethis::use_git_config(user.name = "Andrés Montenegro", user.email = "andresfemole@gmail.com")

# Para crear el git
usethis::use_git()

# Para utilizar github
# Este comando crea automáticamente el repositorio en GIT desde cero
# (hay que generar token y credenciales primero)
usethis::use_github()

# Github solicita autenticación por lo que se debe crear el token
usethis::create_github_token()

# Activar las credenciales de github con el token
gitcreds::gitcreds_set()

# Crear en DESCRIPTION la sección import con las librerias utilizadas
install.packages("attachment")
attachment::att_amend_desc()
