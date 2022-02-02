sidebar_controller = {
    collapsed: false,
    content: "nav",

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
    },

    show: function(content) {
        $(`#${this.content}`).hide()
        $(`#${content}`).show()
        this.content = content
    }
}