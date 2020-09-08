limit <- 1

response_table_old <- response_table
response_table <- response_table_old

# N-P tradeoff
smafit_NP <- sma(trt_P ~ trt_N * lifespan, response_table)

NP_plot <- smafit_NP$groupsummary %>%
  ggplot() +
  geom_abline(aes(slope = Slope,
                  intercept = Int,
                  color = group),
              data = smafit_NP$groupsummary,
              size = 1) +
  geom_point(aes(x = trt_N,
                 y = trt_P,
                 color = lifespan),
             data = response_table,
             alpha = .4,
             size = 1.5) +
  xlim(-limit, limit) +
  ylim(-limit, limit) +
  xlab("Stdzd N Response") +
  ylab("Stdzd P Response") +
  labs(color = "Lifespan") +
  ggtitle("N-P Response Pairs") +
  theme_bw() + 
  theme(legend.position="bottom")+
  theme(plot.title = element_text(size = 12, face = "bold"))

# N-K tradeoff
smafit_NK <- sma(trt_K ~ trt_N * lifespan, response_table)
NK_plot <- smafit_NK$groupsummary %>%
  ggplot() +
  geom_abline(aes(slope = Slope,
                  intercept = Int,
                  color = group),
              data = smafit_NK$groupsummary,
              size = 1) +
  geom_point(aes(x = trt_N,
                 y = trt_K,
                 color = lifespan),
             data = response_table,
             alpha = .4,
             size = 1.5) +
  xlim(-limit, limit) +
  ylim(-limit, limit) +
  xlab("Stdzd N Response") +
  ylab("Stdzd K Response") +
  labs(color = "Lifespan") +
  ggtitle("N-K Response Pairs") +
  theme_bw() + 
  theme(legend.position="bottom")+
  theme(plot.title = element_text(size = 12, face = "bold"))

# P-K tradeoff
smafit_PK <- sma(trt_K ~ trt_P * lifespan, response_table)
PK_plot <- smafit_PK$groupsummary %>%
  ggplot() +
  geom_abline(aes(slope = Slope,
                  intercept = Int,
                  color = group),
              data = smafit_PK$groupsummary,
              size = 1) +
  geom_point(aes(x = trt_P,
                 y = trt_K,
                 color = lifespan),
             data = response_table,
             alpha = .4,
             size = 1.5) +
  xlim(-limit, limit) +
  ylim(-limit, limit) +
  xlab("Stdzd P Response") +
  ylab("Stdzd K Response") +
  labs(color = "Lifespan") +
  ggtitle("P-K Response Pairs")+
  theme_bw() + 
  theme(legend.position="bottom")+
  theme(plot.title = element_text(size = 12, face = "bold"))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(NP_plot)

p3 <- arrangeGrob(arrangeGrob(NP_plot + theme(legend.position="none"),
                              NK_plot + theme(legend.position="none"),
                              PK_plot + theme(legend.position="none"),
                              nrow=1),
                  mylegend, nrow=2,heights=c(10, 1))


ggsave("../../Figures/SMAfits_lifespan.pdf",
       p3,
       height = 4,
       width = 8)
