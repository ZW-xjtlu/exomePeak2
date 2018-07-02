#' @title Dog construction
#'
#' @param personality can be one of the "joyful", "neurotic", "cruel", and "merciful".
#'
#' @description This function defines the constructors for dog class.
Corgi <- function(personality = c("joyful","neurotic",
                                  "cruel","merciful"),
                  weight = 0,
                  length = 0,
                  speed = 0,
                  beauty = 0,
                  satisfaction = 0,
                  fur_amount = 0){

  new("corgi",
      personality = personality[1],
      weight = weight,
      length = length,
      speed = speed,
      beauty = beauty,
      satisfaction = satisfaction,
      fur_amount = fur_amount)

}
