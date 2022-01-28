sidebar_controller = {
    collapsed: false,

    toggle: function() {
        if (this.collapsed) {
            this.collapsed = false
            $("#sidebar-content").show()
            $("#show-hide-button svg").replaceWith(feather.icons["arrow-left"].toSvg())
        } else {
            this.collapsed = true
            $("#sidebar-content").hide()
            $("#show-hide-button svg").replaceWith(feather.icons["arrow-right"].toSvg())
        }
    }
}